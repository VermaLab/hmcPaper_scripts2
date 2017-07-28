SNAME=$1
FQ1=$2
FQ2=$3
REF=/home/kpradha1/projects/hg19/

cat <<PROGRAM_END
#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=40G
#$ -pe smp 17
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu

module load bismark/0.14.5 
module load bowtie2/2.2.3/gcc.4.4.7
module load samtools/1.2/gcc.4.4.7
module load picard-tools/1.92/java.1.8.0_20
