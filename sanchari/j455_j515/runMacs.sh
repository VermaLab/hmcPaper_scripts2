#!/bin/bash
#$ -cwd
#$ -N macs2
#$ -j n
#$ -l h_vmem=9.8G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu

#run mac2 peak finding software on a single target bam
#usually we need a control file, igg or input
#
#usage:
#qsub runMacs.sh bamfile.bam

BAM=$1
NAME=${BAM%.*}
FOL=macs_${NAME}

#module load MACS2/2.1.0/python.2.7.8
module load MACS2/2.1.0-update/python.2.7.8

echo $BAM
echo $NAME
echo $FOL

macs2 callpeak -f BAMPE -t $BAM -n $NAME --outdir $FOL --verbose 3
#macs2 callpeak  -t $BAM -n $NAME --outdir $FOL --call-summits --verbose 3
