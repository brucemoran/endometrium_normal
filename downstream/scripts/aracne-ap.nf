#!/usr/bin/env nextflow

/* run ARACNe-AP in Singularity
*/

params.help = ""

if (params.help) {
  log.info ''
  log.info '----------------------'
  log.info 'NEXTFLOW ARACNe-AP RUN'
  log.info '----------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run aracne-ap.nf \
            --inputCsv "your.gene.matrix.tsv,your.regulators.list,your.seed" \
            --outDir "output/directory/of/choice" \
            -c "aracne-ap.nextflow.config" \
            -with-report "aracne-ap.report.html" \
            -with-timeline "aracne-ap.timeline.html"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --inputCsv     STRING      three objects, full path for first 2: gene expression matrix, tab-delimited, column headers are sample IDs, rownames are gene, with header \'gene\'; gene IDs of regulators (N.B. should match IDs in geneExpMatrix \'gene\' column; third is value to use as seed in threshold calculation'
  log.info '    --outDir     STRING      path to store output'
  log.info ''
  exit 1
}

/* 1.0: calculateThreshold
*/
Channel.fromPath("$params.inputCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [file(row.expmat), file(row.regulators), row.seed] }
       .set { calcthresh }

process calcThresh {

  label 'half_cpu_mem'

  publishDir "$params.outDir", mode: "copy", pattern: "*"

  input:
  set file(expmat), file(regulators), val(seed) from calcthresh

  output:
  set file(expmat), file(regulators), file('*txt') into bootstraping

  script:
  """
  ARACNe-AP ${params.half_javamem} \
            -e $expmat\
            -o ./ \
            -t $regulators \
            -p 1E-8 \
            -s $seed \
            --calculateThreshold
  """
}

/* 2.0: bootstraps
*/
process bootStrap {

  label 'twentieth_cpu_mem'

  publishDir "$params.outDir", mode: "copy", pattern: "*"

  input:
  set file(expmat), file(regulators), file(bootstrap) from bootstraping
  each s from 1..100

  output:
  file('*.txt') into consolidate

  script:
  """
  ARACNe-AP ${params.twentieth_javamem} \
          -e $expmat \
          -o ./ \
          -t $regulators \
          -p 1E-8 \
          -s $s
  """
}

/* 3.0: consolidate bootstraps
*/
process consBoots {

  label 'half_cpu_mem'

  publishDir "$params.outDir", mode: "copy", pattern: "*"

  input:
  file(bootstraps) from consolidate.collect()

  output:
  file('network.txt') into completed

  script:
  """
  ARACNe-AP ${params.half_javamem} \
            -o ./ \
            --consolidate
  """
}
completed.subscribe { println "ARACNe-AP run complete" }
