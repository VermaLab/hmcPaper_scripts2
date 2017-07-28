#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 4
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20

bam=$1
bname="${bam%.*}"

java -jar $(which BuildBamIndex.jar) INPUT=${bam}
java -jar $(which MarkDuplicates.jar) INPUT=${bam} OUTPUT=${bname}_mdup.bam METRICS_FILE=${bname}_mdMetric.txt 
java -jar $(which BuildBamIndex.jar) INPUT=${bname}_mdup.bam


