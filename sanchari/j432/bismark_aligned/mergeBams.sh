#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 4
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
#usage:
#qsub mergeBams.sh outfile.bam infile1.bam infile2.bam infile3.bam
#qsub mergeBams.sh hpne_bs.bam *hpne*bs*sorted*bam
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20

#make a folder for storing intermediate resulst
#otherwise, might run out of space on default /tmp
mkdir -p tmp
chmod 777 tmp

#first argument is the name of the output file
outbam=$1
shift
#the rest of the arguments are the bam files we want to merge

java -jar $(which MergeSamFiles.jar) \
    `echo $@ | sed 's/\(\S\+\)/I=&/g'`\
    ASSUME_SORTED=true\
    TMP_DIR=tmp\
    O=$outbam

