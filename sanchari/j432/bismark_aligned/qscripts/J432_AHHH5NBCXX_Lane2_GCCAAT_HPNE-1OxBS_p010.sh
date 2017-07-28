#!/bin/bash
#$ -cwd
#$ -N bismark_alignment
#$ -j n
#$ -l h_vmem=5.8G
#$ -pe smp 5
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu

module load bismark/0.14.5 
module load bowtie2/2.2.3/gcc.4.4.7
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20

SNAME=J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_p010
FQ1=../trimmed/split/J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS.1_val_1.fq.p010
FQ2=../trimmed/split/J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS.2_val_2.fq.p010
REF=/home/kpradha1/projects/hg19/
bismark --multicore 1 --bowtie2 -p 4 /home/kpradha1/projects/hg19/ -1 ../trimmed/split/J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS.1_val_1.fq.p010 -2 ../trimmed/split/J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS.2_val_2.fq.p010 
java -jar $(which SortSam.jar) INPUT=J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS.1_val_1.fq.p010_bismark_bt2_pe.bam OUTPUT=J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_p010_sorted.bam SO=coordinate CREATE_INDEX=true
echo "done" > J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_p010.status
