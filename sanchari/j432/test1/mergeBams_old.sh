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

mkdir -p tmp
chmod 777 tmp
java -jar $(which MergeSamFiles.jar) \
    `ls *HPNE-1_BS*pe.bam | sed 's/.*/I=&/' | tr '\n' ' '`\
    ASSUME_SORTED=false\
    TMP_DIR=tmp\
    O=test.bam

