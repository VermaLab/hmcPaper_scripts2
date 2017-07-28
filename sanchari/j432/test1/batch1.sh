#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=40G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
module load bismark/0.14.5 
module load bowtie2/2.2.3/gcc.4.4.7
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20

FNAME=test
SNAME=foo
FQ1=${FNAME}_1.fastq.gz
FQ2=${FNAME}_2.fastq.gz
REF=/home/kpradha1/projects/hg19/
#can't use basename with multicore
#bismark --multicore 2 --bowtie2 -p 4 $REF -1 $FQ1 -2 $FQ2  --basename $SNAME
bismark --multicore 4 --bowtie2 -p 4 $REF -1 $FQ1 -2 $FQ2 
#java -jar $(which SortSam.jar) INPUT=${FQ1##*/}_bismark_bt2_pe.bam OUTPUT=${SNAME}_sorted.bam SO=coordinate
java -jar $(which SortSam.jar) INPUT=${FQ1}_bismark_bt2_pe.bam OUTPUT=${SNAME}_sorted.bam SO=coordinate
java -jar $(which BuildBamIndex.jar) INPUT=${SNAME}_sorted.bam
java -jar $(which MarkDuplicates.jar) INPUT=${SNAME}_sorted.bam OUTPUT=${SNAME}_sorted_md.bam METRICS_FILE=${SNAME}_mdMetric.txt
java -jar $(which BuildBamIndex.jar) INPUT=${SNAME}_sorted_md.bam


