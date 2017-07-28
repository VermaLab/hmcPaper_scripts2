#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
module load samtools/1.2/gcc.4.4.7
module load bamUtil/1.0.13/gcc.4.4.7

bam=$1
bname="${bam%.*}"

bam clipOverlap --in $bam --out ${bname}_clip.bam
samtools index ${bname}_clip.bam


