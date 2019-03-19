#!/usr/bin/env nextflow

/* QC, trimming of paired-end Illumina data;
* fastp QC; bbduk trim by quality into single fastq
* align with STAR+markDuplicates; psueodalign with kallisto
* metrics from picard suite, visualised by multiQC; combine featureCounts into matrix
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '--------------------------------------------'
  log.info 'NEXTFLOW RNASEQ QC, TRIM, ALIGN, PSEUDOALIGN'
  log.info '--------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run hg38.RNAseq_STAR-kallisto.simg.nf \
            --sampleCsv "/full/path/to/readsDir/sampleInput.csv" \
            --projectBase "/full/path/to/base" \
            --fa "/path/to/reference/your.fa" \
            --gtf "/path/to/reference/your.gtf" \
            --STAR "/path/to/reference/STAR/your" \
            --kallisto "/path/to/reference/your.kallisto" \
            --rRNA "/path/to/reference/your.dict.interval_list" \
            --refFlat "/path/to/reference/your.gtf.refFlat" \
            --noStar \
            -c "RNAseq_STAR-kallisto.simg.nextflow.config" \
            -with-report "RNAseq_STAR-kallisto.report.html" \
            -with-timeline "RNAseq_STAR-kallisto.timeline.html"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --sampleCsv     STRING      /path/to/csv file with header: sampleID,stranded,reads1,reads2; NB no read2 required, omission runs sample as single-ended'
  log.info '    --projectBase     STRING      base directory for project in which runID is made'
  log.info '    --fa     STRING      human fasta (faidx\'d)'
  log.info '    --gtf    STRING      human gtf'
  log.info '    --STAR    STRING      human STAR index dir'
  log.info '    --noSTAR    STRING      simple flag to not run STAR or any of the concomitant processes'
  log.info '    --kallisto     STRING      human kallisto index'
  log.info '    --rRNA     STRING      human rRNA for picard-tools QC'
  log.info ''
  exit 1
}

params.outDir = "$params.projectBase/analysis/RNAseq"

//for PE, SE: https://groups.google.com/forum/#!topic/nextflow/_ygESaTlCXg
Channel.fromPath("$params.sampleCsv", type: 'file')
    .splitCsv( header: true )
    .map { row -> [row.sampleID, row.stranded, [file(row.read1), file(row.read2)] - null ] }
    .set { bbduking }

/* 1.0: Input trimming
 * NB each fastq in data/fastq has own dir, name is mirrored in analysis/
 * here we start the data throughput
 */
process bbduk {

  label 'eighth_cpu_mem'
  publishDir "$params.outDir/$sampleID/bbduk", mode: "copy", pattern: "*.bbduk.runstats.txt"

  input:
  set val(sampleID), val(stranded), file(reads) from bbduking

  output:
  set val(sampleID), file(reads) into fastping
  set val(sampleID), val(stranded), file('*bbduk.*.fastq.gz') into (staring, kallistoing)

  script:
  def single = reads instanceof Path
  if( !single ) {
    """
    READ1=\$(ls *.R1.*); READ2=\$(ls *.R2.*)
    { bbduk.sh ${params.quarter_javamem} \
        in1=\$READ1 out1=$sampleID".bbduk.R1.fastq.gz" \
        in2=\$READ2 out2=$sampleID".bbduk.R2.fastq.gz" \
        k=${params.bbdkmerx} mink=${params.bbdmink} \
        hdist=1 ktrim=r trimq=${params.bbdqtrim} \
        qtrim=rl maq=20 ref=${params.bbmapAdapters} tpe tbo \
        stats=$sampleID".bbduk.adapterstats.txt" overwrite=T

    } 2>&1 | tee $sampleID".bbduk.runstats.txt"
    """
  }
  else {
    """
    { bbduk.sh ${params.quarter_javamem} \
        in1=$reads out1=$sampleID".bbduk.R1.fastq.gz" \
        k=${params.bbdkmerx} mink=${params.bbdmink} \
        hdist=1 ktrim=r trimq=${params.bbdqtrim} \
        qtrim=rl maq=20 ref=${params.bbmapAdapters} tpe tbo \
        stats=$sampleID".bbduk.adapterstats.txt" overwrite=T

    } 2>&1 | tee $sampleID".bbduk.runstats.txt"
    """
  }
}


/* 1.1: Input QC
 */
process fastp {

  label 'eighth_cpu_mem'

  input:
  set val(sampleID), file(reads) from fastping

  output:
  file('*') into completedfastp
  file('*.json') into fastp_multiqc

  script:
  def single = reads instanceof Path
  if( !single ) {
    """
    READ1=\$(ls *.R1.*); READ2=\$(ls *.R2.*)
    fastp -w ${task.cpus} -j $sampleID".fastp.json" --in1 \$READ1 --in2 \$READ2
    """
  }
  else {
    """
    fastp -w ${task.cpus} -j $sampleID".fastp.json" --in1 $reads
    """
  }
}

/* 2.0: STAR align and count hg38
*/
process star {

  label 'full_cpu_mem'
  publishDir "$params.outDir/$sampleID/STAR", mode: "copy", pattern: "${sampleID}.*[!bam,!bai]"

  input:
  set val(sampleID), val(stranded), file(reads) from staring

  output:
  file('*') into completedstar
  set val(sampleID), file('*.Aligned.sortedByCoord.out.bam') into dup_marking
  file("${sampleID}.Log.final.out") into starLOG_multiqc

  when:
  !params.noStar

  script:
  """
  BAMsortRAM=\$(echo ${task.memory} | sed 's/ G/000000000/' | \
    perl -ane '\$o=\$F[0]-1000000000; print "\$o\\n";')
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE=\$(echo ID:${sampleID}\\" \\"PL:ILLUMINA\\" \\"SM:${sampleID}\\" \\"DS:RNAseq\\" \\"CN:UNK\\" \\"LB:LANE_X\\" \\"DT:\$DATE)

  STAR --runThreadN ${task.cpus} \
     --genomeDir ${params.STAR} \
     --readFilesIn $reads \
     --readFilesCommand "zcat" \
     --outFileNamePrefix $sampleID"." \
     --outSAMmode Full \
     --outSAMstrandField intronMotif \
     --outSAMattributes All \
     --outFilterType BySJout \
     --seedSearchStartLmaxOverLread 0.5 \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode TranscriptomeSAM GeneCounts \
     --twopassMode Basic \
     --limitBAMsortRAM \$BAMsortRAM \
     --outSAMattrRGline \$RGLINE
  """
}

/* 2.1: Duplicate Marking
*/
process mrkdup {

  label 'full_cpu_mem'
  publishDir "$params.outDir/$sampleID/STAR", mode: "copy", pattern: "*"

  input:
  set sampleID, file(bam) from dup_marking

  output:
  set val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into (featurecounting, picardmetricing)
  file('*.md.metrics.txt') into (completedmarkduplicates, markdup_multiqc)

  when:
  !params.noStar

  script:
  """
  #! /bin/bash

  samtools index $bam

  {

  OUTBAM=\$(echo $bam | sed 's/bam/md.bam/')
  OUTMET=\$(echo $bam | sed 's/bam/md.metrics.txt/')

  picard-tools ${params.full_javamem} \
    MarkDuplicates \
    TMP_DIR=./ \
    INPUT=$bam \
    OUTPUT=/dev/stdout \
    COMPRESSION_LEVEL=0 \
    METRICS_FILE=\$OUTMET \
    REMOVE_DUPLICATES=FALSE \
    ASSUME_SORTED=TRUE \
    VALIDATION_STRINGENCY=LENIENT | samtools view -Shb - > \$OUTBAM

  samtools index \$OUTBAM > \$OUTBAM".bai"

  } 2>&1 | tee $sampleID".md.log.txt"
  """
}

/* 2.2: FeatureCount
*/
FESSES = Channel.from( 0, 1, 2 )

process featco {

  label 'quarter_cpu_mem'
  publishDir "$params.outDir/$sampleID/featurecounts", mode: "copy", pattern: "*"

  input:
  set sampleID, file(bam), file(bai) from featurecounting
  each fess from FESSES

  output:
  file('*') into completedfeaturecount
  file("${sampleID}.Aligned.sortedByCoord.out.md.featco.*.counts") into featurecount_combine
  file("${sampleID}.Aligned.sortedByCoord.out.md.featco.*.tab.summary") into featurecount_multiqc

  when:
  !params.noStar

  script:
  """
  #! /bin/bash

  OUTPUT=\$(echo $bam | sed "s/bam/featco.s${fess}.tab/")
  OUTCOU=\$(echo $bam | sed "s/bam/featco.s${fess}.counts/")

  {
  featureCounts \
    -F "SAF" \
    -a "saf.saf" \
    -o \$OUTPUT \
    -s ${fess} \
    -T 10 \
    --primary \
    --ignoreDup \
    -p -P -D 100000000 \
    $bam

  cat \$OUTPUT | perl -ane 'if((\$F[0]=~m/^Geneid/) || (\$F[0]=~m/^#/)){next;}else{chomp;print "\$F[0]\\t\$F[6]\\n"}' >> \$OUTCOU
  } 2>&1 | tee > $sampleID".featurecounts.s${fess}.log.txt";
  """
}

/* 2.3: Picard Metricst
*/
process picmet {

  label 'quarter_cpu_mem'
  publishDir "$params.outDir/$sampleID/picard", mode: "copy", pattern: "${sampleID}.*"

  input:
  set sampleID, file(bam), file(bai) from picardmetricing

  output:
  file('*') into completedpicardmetric
  file('*.CollectRnaSeqMetrics.*') into picmet_rnaseq_multiqc
  file('*.AlignmentSummaryMetrics.*') into picmet_alignsummary_multiqc
  file('*.CollectMultipleMetrics.*') into picmet_collectmultiple_multiqc
  file('*.est_lib_complex_metrics.*') into picmet_estlibcomplex_multiqc
  file('*.insert_size_metrics.*') into picmet_insertsize_multiqc

  when:
  !params.noStar

  script:
  """
  #! /bin/bash
  {
  ##make rRNA.interval_list as BAM
  grep -v "@" ${params.rRNA} | perl -ane '@s=split(//,\$F[0]);if(scalar(@s)>5){next;}else{print \$_;}' > rRNA.bed
  picard-tools BedToIntervalList \
      I=rRNA.bed \
      O=rRNA.interval_list \
      SD=$bam

  picard-tools CollectRnaSeqMetrics \
    TMP_DIR=./ \
    REF_FLAT=${params.refFlat} \
    RIBOSOMAL_INTERVALS=rRNA.interval_list \
    STRAND_SPECIFICITY="NONE" \
    CHART_OUTPUT=$bam".norm-pos_vs_cov.pdf" \
    INPUT=$bam \
    OUTPUT=$bam".CollectRnaSeqMetrics.txt"

    picard-tools CollectAlignmentSummaryMetrics \
    I=$bam \
    O=$sampleID".AlignmentSummaryMetrics.txt" \
    TMP_DIR=./ \
    R=${params.fa}

  picard-tools CollectMultipleMetrics \
    I=$bam \
    O=$sampleID".CollectMultipleMetrics.txt" \
    TMP_DIR=./ \
    R=${params.fa}

  picard-tools CollectSequencingArtifactMetrics \
    I=$bam \
    O=$sampleID".artifact_metrics.txt" \
    TMP_DIR=./ \
    R=${params.fa}

  picard-tools EstimateLibraryComplexity \
    I=$bam \
    O=$sampleID".est_lib_complex_metrics.txt" \
    TMP_DIR=./

  picard-tools CollectInsertSizeMetrics \
    I=$bam \
    O=$sampleID".insert_size_metrics.txt" \
    H=$bam".histogram.pdf" \
    TMP_DIR=./

  } 2>&1 | tee > $sampleID".picard.metrics.log.txt"
  """
}


/* 2.4: Kallisto
*/
process kallisto {

  label 'eighth_cpu_mem'
  publishDir "$params.outDir/$sampleID/kallisto", mode: "copy", pattern: "*"

  input:
  set val(sampleID), val(stranded), file(reads) from kallistoing

  output:
  file('*') into completedkallisto
  file('*.kallisto.log.txt') into kallisto_multiqc

  script:
  def single = reads instanceof Path
  if( !single ) {
    """
      { kallisto quant \
          -t 10 -b 100 $stranded -i ${params.kallisto} -o ./ $reads
      } 2>&1 | tee > $sampleID".kallisto.log.txt"
    """
  }
  else {
    """
      { kallisto quant \
          --single -l 200 -s 30 -t 10 -b 100 $stranded -i ${params.kallisto} -o ./ $reads
      } 2>&1 | tee > $sampleID".kallisto.log.txt"
    """
  }
}

/* 3.0: MultiQC
 * NB bbduk_multiqc output not available for multic yet 130718
*/
fastp_multiqc.mix(kallisto_multiqc).set { first_multiqc }

if( params.noStar == true ){
  starLOG_multiqc.mix( featurecount_multiqc, markdup_multiqc, picmet_rnaseq_multiqc, picmet_alignsummary_multiqc, picmet_collectmultiple_multiqc, picmet_estlibcomplex_multiqc, picmet_insertsize_multiqc ).set { second_multiqc }
}

first_multiqc.mix(second_multiqc).set { multiqc_multiqc }

process multiqc {

  label 'eighth_cpu_mem'
  publishDir "$params.outDir/multiqc", mode: "copy", pattern: "*"

  input:
  file ('*') from multiqc_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  multiqc -n multiqc_report.html -c ${params.multiqcconfig} -f .
  """
}

/* 3.1: Combine FeatureCounts
*/
ESSES = Channel.from( 0, 1, 2 )

process cmbfco {

  label 'quarter_cpu_mem'
  publishDir "$params.outDir/DESeq2", mode: "copy", pattern: "*"

  input:
  file(counts) from featurecount_combine.collect()
  each ess from ESSES

  output:
  file('featureCounts.s*.tsv') into completedfeaturecountscombine

  when:
  !params.noStar

  script:
  """
  #! /bin/bash
  FRST=\$(ls *s${ess}.counts | sort -V | head -n1)
  cat \$FRST > 1
  ls *s${ess}.counts | sort -V | tail -n+2 | \
    while read FILE; do
      join 1 \$FILE > 2
      mv 2 1
    done
  echo -n "ensembl_gene_id" > 2
  ls *s${ess}.counts | sort -V | while read HEADER; do
    HDR=\$(echo \$HEADER | cut -d "." -f1)
    echo -ne "\\t"\$HDR >> 2
  done
  echo "" >> 2
  cat 1 >> 2
  sed 's/\\s */\\t/g' 2 > featureCounts.s${ess}.tsv;
  """
}
