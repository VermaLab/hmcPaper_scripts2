#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=15G
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
module load bismark/0.14.5 
module load bowtie2/2.2.3/gcc.4.4.7
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20

FNAME=J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS
SNAME=hpne_bs
FQ1=${FNAME}.1_val_1.fq.gz
FQ2=${FNAME}.2_val_2.fq.gz
REF=/home/kpradha1/projects/hg19/
bismark --bowtie2 $REF -1 $FQ1 -2 $FQ2  --basename $SNAME
java -jar $(which SortSam.jar) INPUT=${SNAME}_pe.bam OUTPUT=${SNAME}_sorted.bam SO=coordinate
java -jar $(which BuildBamIndex.jar) INPUT=${SNAME}_sorted.bam
java -jar $(which MarkDuplicates.jar) INPUT=${SNAME}_sorted.bam OUTPUT=${SNAME}_sorted_md.bam METRICS_FILE=${SNAME}_mdMetric.txt
java -jar $(which BuildBamIndex.jar) INPUT=${SNAME}_sorted_md.bam


