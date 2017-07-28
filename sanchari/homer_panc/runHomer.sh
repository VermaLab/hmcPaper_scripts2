#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 4
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu

#add thep path to homer and seqlogo
#PATH=$PATH:/home/kpradha1/programs/homer/bin
PATH=$PATH:/home/kpradha1/programs/weblogo
module load blatSuite/34
module load HOMER/4.7

outname1=$1
peaks1=$2


findMotifsGenome.pl $peaks1 hg19 $outname1 -len 0 -preparsedDir /home/kpradha1/programs/homer_preparse -preparse -p 8

