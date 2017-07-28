#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=15G
#$ -S /bin/bash
#$ -m e
#$ -M kith.pradhan@einstein.yu.edu
#
# do_bismark.sh
# Brent Calder <brent.calder@einstein.yu.edu>
#
# Convienience script to take fastq.gz input, run trim_galore,
# fastqc, split input, run Bismark, then merge, sort, and index
# output into BAM.
#


set -o errexit -o nounset -o allexport -o pipefail

# user options

GENOME_FOLDER=/home/kpradha1/projects/hg19/
GENOME=mm9.fa # not used, this is the original 
# genome name at ../fasta, relative to GENOME_FOLDER

TRIM_GALORE_OPTS="" # mainly for new adapter sequence or clips
FASTQC_OPTS="" # optional FASTQC options
BISMARK_OPTS="" # see documentation, 'non-directional' set below (DIRECTIONS) 
BOWTIE2_OPTS="" # bowtie2 options for bismark

# number of reads per split work file
READS_PER_FASTQ=4000000
# DIRECTIONS must be 2 or 4 (eg. 4 = "--non_directional")
DIRECTIONS=4
# NTHREADS is the number of threads per DIRECTION
NTHREADS=2
# DIRECTIONS*THREADS = num_proc

########################
# DO NOT EDIT BELOW HERE

USAGE="Usage:\t\tqsub do_bismark.sh forward_read.fastq.gz [reverse_read.fastq.gz]\n"
ADMON="reads must be in format: name.fastq.gz\n"
PAIRED=0
REVERSE_READ=""
MYID=$$-${RANDOM}

if [ "$DIRECTIONS" -eq "4" ]; then
    BISMARK_OPTS="${BISMARK_OPTS} --non_directional"
else
    # default back to directional
    DIRECTIONS=2
fi

if [ "$#" -eq "2" ]; then
    PAIRED=1
    REVERSE_READ=$2
    REVERSE=`basename ${REVERSE_READ/.fastq.gz/}`
    TRIM_GALORE_OPTS="${TRIM_GALORE_OPTS} --paired"
    if [ "${#REVERSE_READ}" -eq "${#REVERSE}" ]; then
        USAGE="${USAGE}${ADMON}"
        echo -e "$USAGE" >&2
        exit 2
    fi
else
    if [ "$#" -ne "1" ]; then
        echo -e "$USAGE" >&2
        exit 1
    fi
fi

FORWARD_READ=$1
FORWARD=`basename ${FORWARD_READ/.fastq.gz/}`

if [ "${#FORWARD_READ}" -eq "${#FORWARD}" ]; then
    USAGE="${USAGE}${ADMON}"
    echo -e "$USAGE" >&2
    exit 2
fi

##################
# REQUIRED MODULES
##################
#source /apps1/modules/init/bash

module load trim_galore/0.3.7/gcc.4.4.7
module load bismark/0.14.5
module load FastQC/0.11.4/java.1.8.0_20
module load picard-tools/1.92/java.1.8.0_20

# might be necessary to set up qsub
if [ -e /etc/profile.d/sge-binaries.sh ]; then
    source /etc/profile.d/sge-binaries.sh
fi

# output dirs
mkdir -p fastqc split bismark trimmed tmp

# trim illumina adapters and -q 20 by default
trim_galore --dont_gzip -o trimmed --fastqc --fastqc_args "${FASTQC_OPTS} -o ./fastqc" $TRIM_GALORE_OPTS $FORWARD_READ $REVERSE_READ

# split files based on READS_PER_FASTQ
if [ "$PAIRED" -eq "0" ]; then
    split -a 3 -l $((${READS_PER_FASTQ}*4)) trimmed/${FORWARD}_trimmed.fq split/${FORWARD}.split.
else
    if [ "$PAIRED" -eq "1" ]; then
        split -a 3 -l $((${READS_PER_FASTQ}*4)) trimmed/${FORWARD}_val_1.fq split/${FORWARD}.split.
        split -a 3 -l $((${READS_PER_FASTQ}*4)) trimmed/${REVERSE}_val_2.fq split/${REVERSE}.split.
    fi
fi

# run bismark on the split files
for x in `ls -1 split/${FORWARD}.split.*`; do
    READSTR=${x}
    suf=${x/*split./}
    if [ "$PAIRED" -eq "1" ]; then
        second=split/${REVERSE}.split.${suf}
        READSTR="-1 ${x} -2 ${second}"
    fi
    qsub -cwd -l p=$((${NTHREADS}*${DIRECTIONS})) -l mf=30g -V -N bismark-${MYID} -j y -b y bismark --temp_dir tmp -o bismark --bowtie2 -p $((${NTHREADS}*${DIRECTIONS})) $BISMARK_OPTS $BOWTIE2_OPTS $GENOME_FOLDER $READSTR
done

# use picard to merge, sort, create bam and index split SAM upon completion.
SAMS=""
if [ "$PAIRED" -eq "0" ]; then
    SAMS=`ls -1 split/${FORWARD}.split.* | sed 's|^split|I=bismark|g' | sed 's/$/_bismark_bt2.sam/g' | tr "\n" " "`
else
    SAMS=`ls -1 split/${FORWARD}.split.* | sed 's|^split|I=bismark|g' | sed 's/$/_bismark_bt2_pe.sam/g' | tr "\n" " "`
fi
BAM=${FORWARD}.bam

qsub -S /bin/bash -hold_jid bismark-${MYID} -cwd -V -l mf=4g -j y -N merge -b y java -jar $PICARD_ROOT/MergeSamFiles.jar $SAMS O=$BAM AS=false SO=coordinate TMP_DIR=./tmp CREATE_INDEX=true

