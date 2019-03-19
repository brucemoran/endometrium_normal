#!/usr/bin/env nextflow

/* QC, trimming of paired-end Illumina data;
* fastqc; bbduk trim by quality into single fastq
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
  log.info 'nextflow run RNAseq_STAR-kallisto_references.simg.nf \
            --refBase "/full/path/to/base" \
            --fa "ftp://ftp.ensembl.org/pub/release-94/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz" \
            --gtf "ftp://ftp.ensembl.org/pub/release-94/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.94.gtf.gz" \
            --cdna "ftp://ftp.ensembl.org/pub/release-94/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz" \
            --vcf "ftp://ftp.ensembl.org/pub/release-95/variation/vcf/rattus_norvegicus/rattus_norvegicus.vcf.gz" \
            --set "rat" \
            -with-report "RNAseq_STAR-kallisto_references.report.html" \
            -with-timeline "RNAseq_STAR-kallisto_references.timeline.html" \
            -c RNAseq_STAR-kallisto_references.simg.nextflow.config'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --refBase     STRING      base directory for reference data'
  log.info '    --fa     STRING      direct link to fasta'
  log.info '    --gtf     STRING      direct link to gtf'
  log.info '    --cdna     STRING      direct link to cdna'
  log.info '    --vcf     STRING      direct link to vcf for use in GATK4'
  log.info '    --set     STRING      over-ride above definitions, using a set within config, named as per here'
  log.info ''
  exit 1
}

params.outDir = "$params.refBase/RNAseq_STAR-kallisto_references"

/*-1.0 test for set and use accordingly
*/
if ("$params.set" != "false"){
  params.fa = "\$${params.set}fa"
  params.gtf = "\$${params.set}gtf"
  params.cdna = "\$${params.set}cdna"
  params.vcf = "\$${params.set}vcf"
}

/* 0.0 inputs
*/
process unzpfa {

  publishDir "$params.outDir", mode: "copy", pattern: "*a"

  output:
  file('*a') into (star_fa, rrna_fa, gatk4_fa)

  script:
  """
  wget ${params.fa}
  gunzip *gz
  """
}

process unzpgtf {

  publishDir "$params.outDir", mode: "copy", pattern: "*gtf"

  output:
  file('*gtf') into (star_gtf, refflat_gtf, rrna_gtf, gatk4_gtf, saf_gtf)

  script:
  """
  wget ${params.gtf}
  gunzip *gz
  """
}

process unzpcdna {

  output:
  file('*') into (kallisto_cdna)

  script:
  """
  wget ${params.cdna}
  """
}

process unzpvcf {

  publishDir "$params.outDir", mode: "copy", pattern: "*"

  output:
  file('*') into (gatk4_vcf)

  script:
  """
  wget ${params.vcf}
  BASENAME=\$(echo \$(basename ${params.vcf}) | sed 's/\\.gz//')
  gunzip \$(basename ${params.vcf})
  gatk IndexFeatureFile -F \$BASENAME
  """
}

Channel.from(74, 99)
        .set { star_sjdb }

/* 1.0: STAR geneomeGenerate
 */
process stargg {

  label 'full_cpu_mem'
  publishDir "$params.outDir", mode: "copy", pattern: "*"

  input:
  file(fa) from star_fa
  file(gtf) from star_gtf
  val(sjdb) from star_sjdb

  output:
  file('*') into completedstargg

  script:
  """
  #! /bin/bash
  mkdir -p STAR_$sjdb

  RAM=\$(echo ${task.memory} | sed 's/ G/000000000/')
  STAR --runMode genomeGenerate \
    --genomeDir STAR_$sjdb \
    --genomeFastaFiles $fa \
    --sjdbGTFfile $gtf \
    --sjdbOverhang $sjdb \
    --runThreadN ${task.cpus} \
    --limitGenomeGenerateRAM \$RAM
  """
}

/* 2.0: refFlat conversion of GTF
 */
process reflat {

  label 'eighth_cpu_mem'
  publishDir "$params.outDir", mode: "copy", pattern: "refFlat*"

  input:
  file(gtf) from refflat_gtf

  output:
  file('*') into completedrefflat

  script:
  """
  #! /bin/bash
  wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
  chmod a+x ./gtfToGenePred
  mkdir -p refFlat
  cd refFlat
  ../gtfToGenePred -includeVersion ../$gtf $gtf".refFlat"
  perl -ane 'print "\$F[0]\\t\$_";' $gtf".refFlat" > 1
  mv 1 $gtf".refFlat"
  """
}

/* 3.0: rRNA_intervalList
 */
process rrna {

  label 'eighth_cpu_mem'
  publishDir "$params.outDir", mode: "copy", pattern: "rRNA*"

  input:
  file(fa) from rrna_fa
  file(gtf) from rrna_gtf

  output:
  file('*') into completedrRNA

  script:
  """
  #! /bin/bash

  mkdir -p rRNA
  cd rRNA
  picard-tools CreateSequenceDictionary \
    R=../$fa \
    O=../$fa".dict"

  ##https://www.biostars.org/p/67079/
  cat ../$fa".dict" > $fa".dict.rRNA.interval_list"

  grep 'gene_biotype "rRNA"' ../$gtf | \
    awk '\$3 == "transcript"' | \
    cut -f1,4,5,7,9 | \
    perl -lane '/transcript_id "([^"]+)"/ or die "no transcript_id on \$.";
      print join "\\t", (@F[0,1,2,3], \$1);' | \
    sort -k1V -k2n -k3n >> $fa".dict.rRNA.interval_list"
  """
}

/* 4.0: kallisto index
 */
process kalsto {

  label 'full_cpu_mem'
  publishDir "$params.outDir", mode: "copy", pattern: "kallisto"

  input:
  file(cdna) from kallisto_cdna

  output:
  file('*') into completedkal

  script:
  """
  #! /bin/bash

  mkdir -p kallisto
  cd kallisto
  kallisto index -i $cdna".kallisto" ../$cdna
  """
}

/*5 GATK4 SNV calling for tumour-only required files
*/
process gatk4snv {

  label 'full_cpu_mem'
  publishDir path: "$params.outDir", mode: 'copy'

  input:
  file(gtf) from gatk4_gtf
  file(fa) from gatk4_fa

  output:
  file('*') into gatk4all

  """
  #! /bin/bash

  ##faidx fasta
  samtools faidx $fa

  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
  REFD=\$(echo $fa | sed -r 's/(.fa).*/.dict/')
  picard-tools CreateSequenceDictionary REFERENCE=$fa OUTPUT=\$REFD

  perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' \$REFD > seq.dict.chr-size.dict

  ##make GTF bed for 'gene'
  GENEINTLIST=\$(echo $gtf | sed 's/\\.gtf/\\.genes\\.interval_list/')
  cat $gtf | perl -ane 'if(\$F[2] eq "gene"){print "\$F[0]\\t\$F[3]\\t\$F[4]\\n";}' | \
    perl -ane 'if((\$F[0]=~m/^chr/) || (\$F[0]=~m/^[0-9]/)){\$_=~s/chr//g;print \$_;}' > genes.bed

  ##test for 'chr' in fasta dictionary
  ##always make interval list so we are in line with fasta
  CHRT=\$(cat seq.dict.chr-size.dict | head -n1 | \
    perl -ane '@s=split(/:/,\$F[1]);if(\$s[1]=~m/^chr/){print "CHR\\n"};')

  if [[ \$CHRT == "" ]];then
    cat seq.dict.chr-size.dict | perl -ane 'print "\$F[0]\\n";' | \
    while read CHR; do
      export CHR;
      sort -V genes.bed | perl -ane 'if(\$F[0] eq \$ENV{CHR}){print \$_;}else{next;}'
    done > 1
    picard-tools BedToIntervalList I=1 O=\$GENEINTLIST SD=\$REFD
    rm 1
  else
    picard-tools BedToIntervalList I=genes.bed O=\$GENEINTLIST SD=\$REFD
  fi

  ##tabix
  bgzip \$GENEINTLIST
  gunzip -c \$GENEINTLIST".gz" > \$GENEINTLIST
  tabix \$GENEINTLIST".gz"
  """
}
gatk4all.subscribe { println "Made: " + it }

/*6 SAF for featureCounts
*/
process gtfSaf {

  publishDir path: "$params.outDir", mode: 'copy'

  input:
  file(gtf) from saf_gtf

  output:
  file('*.saf') into saf_out

  """
  #! /bin/bash

  ##make SAF format from gtf
  SAF=\$(echo $gtf | sed 's/\\.gtf/\\.saf/')
  echo -e "GeneID\\tChr\\tStart\\tEnd\\tStrand" > \$SAF
  grep -v "#" $gtf | perl -ane 'chomp;if(\$F[2] eq "gene"){\$g=\$F[9];\$g=~s/[";]//g;print "\$g\\t\$F[0]\\t\$F[3]\\t\$F[4]\\t\$F[6]\\n";}' >> \$SAF
  """
}
saf_out.subscribe { println "Made: " + it }
