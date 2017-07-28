#!/bin/bash
#constructs the qsub script to process a chunk of reads
#writes the script to stdout
#pipe it into qsub if you want to submit it
#
#usage:
#makeBismarkBatch.sh hpneBs reads1.fastq.gz reads2.fastq.gz | qsub

SNAME=$1
FQ1=$2
FQ2=$3
REF=/home/kpradha1/projects/hg19/

cat <<PROGRAM_END
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

SNAME=$1
FQ1=$2
FQ2=$3
REF=/home/kpradha1/projects/hg19/
bismark --multicore 1 --bowtie2 -p 4 $REF -1 $FQ1 -2 $FQ2 
java -jar \$(which SortSam.jar) INPUT=${FQ1##*/}_bismark_bt2_pe.bam OUTPUT=${SNAME}_sorted.bam SO=coordinate CREATE_INDEX=true
echo "done" > $SNAME.status
PROGRAM_END

