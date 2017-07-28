#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
module load samtools/1.2/gcc.4.4.7
module load python/2.7.8/gcc.4.4.7

REF=/home/kpradha1/projects/hg19/hg19.fa
bam=$1
bname="${bam%.*}"


samtools mpileup -d1000000 -Q0 -B -f $REF $bam | \
  ../mpileup2methylation.py -i - >| ${bname}.bedGraph



