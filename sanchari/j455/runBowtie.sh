#!/bin/bash
#$ -cwd
#$ -N bowtie_alignment
#$ -j n
#$ -l h_vmem=5.8G
#$ -pe smp 2
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu

#process a pair of fastqs
#
#usage:
#qsub runBowtie.sh aligned_name reads1.fastq.gz reads2.fastq.gz

SNAME=$1
FQ1=$2
FQ2=$3
REF=../../hg19/hg19

module load bowtie2/2.2.3/gcc.4.4.7
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20


bowtie2 -t -p 2 -q --local -x $REF -1 <(zcat $FQ1) -2 <(zcat $FQ2) 2>| summary_${SNAME}.txt | samtools view -bS - >| ${SNAME}.bam
java -jar $(which SortSam.jar) INPUT=${SNAME}.bam OUTPUT=${SNAME}_sorted.bam SO=coordinate CREATE_INDEX=true
echo "done" >| $SNAME.status

