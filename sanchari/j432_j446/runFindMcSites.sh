#!/bin/bash
#$ -cwd
#$ -N findMcSites
#$ -j n
#$ -l h_vmem=1G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu

#run the streaming r script that finds mc sites
#
#usage:
#qsub runFindMcSites.sh HPNE_BS_1.txt.gz HPNE_OxBS_1.txt.gz mcSites/stream_mcSites_HPNE_1.txt

BSFILE=$1
OXBSFILE=$2
OUTFILE=$3

module load R/3.3.0/gcc.4.4.7

echo $BSFILE
echo $OXBSFILE
echo $OUTFILE

./findMcSites_stream.r $BSFILE $OXBSFILE $OUTFILE
