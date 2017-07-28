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

BAM1=$1
BAM2=$2
NAME=C${BAM1%.*}_T${BAM2%.*}
FOL=macs_${NAME}

#module load MACS2/2.1.0/python.2.7.8
module load MACS2/2.1.0-update/python.2.7.8

echo $BAM1
echo $BAM2
echo $NAME
echo $FOL

macs2 callpeak -f BAMPE -c $BAM1 -t $BAM2 -n $NAME --outdir $FOL --verbose 3
