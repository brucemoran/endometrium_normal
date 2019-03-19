#! /bin/bash

## for running on genome.ucd.ie

BASEDIR=$1
TAG=$2
DATADIR=$BASEDIR"/ARACNe-AP/"$TAG

if [[ ! -d $DATADIR/bootstraps/qsub ]];then
  mkdir -p $DATADIR/bootstraps/qsub
fi

echo "
## load required (NB ARANCe loads Java8)
module load 'ARACNe-AP/0.0.1'
module load 'R/3.4.0'

## single threshold
ARACNe-AP -Xmx80G \
          -e $DATADIR/dataset.vst.tsv \
          -o $DATADIR/bootstraps\
          -t $DATADIR/signature.genes.txt \
          -p 1E-8 \
          -s 123456789 \
          --calculateThreshold

## networks for 100 bootstraps
for x in {1..100};do
  echo -e \"#! /bin/bash
#PBS -N ARCNe_\${x}
#PBS -q batch
#PBS -l nodes=1:ppn=10
module load 'ARACNe-AP/0.0.1'
module load 'R/3.4.0'

ARACNe-AP -Xmx24G \
          -e $DATADIR/dataset.vst.tsv \
          -o $DATADIR/bootstraps \
          -t $DATADIR/signature.genes.txt \
          -p 1E-8 \
          -s \$x
  \" > $DATADIR/bootstraps/qsub/bootstrap_\${x}.qsub.sh
  qsub $DATADIR/bootstraps/qsub/bootstrap_\${x}.qsub.sh
done 2>&1 | tee > $DATADIR/bootstraps/qsub/HOLDJOBIDS

##sleep, then get all ARACNe processes queued, tell this last consolidate not to run until they complete
sleep 10

echo -e \"#! /bin/bash
#PBS -N ARCNe_C
#PBS -q batch
#PBS -l nodes=1:ppn=40
module load 'ARACNe-AP/0.0.1'
module load 'R/3.4.0'

ARACNe-AP -Xmx80G \
          -o $DATADIR/bootstraps \
          --consolidate

mv $DATADIR/bootstraps/network.txt $DATADIR/bootstraps/network_ready.txt
\" > $DATADIR/bootstraps/qsub/consolidate.qsub.sh

HOLDJIDS=\$(cat $DATADIR/bootstraps/qsub/HOLDJOBIDS | perl -ane 'chomp;if(\$. == 1){print \$F[0];}else{print ":\$F[0]";}')
qsub -W depend:afterok:$HOLDJIDS $DATADIR/bootstraps/qsub/consolidate.qsub.sh
" > $DATADIR/bootstraps/$TAG.ARACNe_AP_run.sh
