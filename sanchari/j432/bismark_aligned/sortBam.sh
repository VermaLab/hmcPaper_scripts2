#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
module load samtools/1.2/gcc.4.4.7

bam=$1
bname="${bam%.*}"

samtools sort $bam ${bname}_sorted
samtools index ${bname}_sorted.bam


