#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=4G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#   #$ -M kith.pradhan@einstein.yu.edu
module load python/2.7.8/gcc.4.4.7
module load numpy/1.9.0/python.2.7.8

inFile=$1
outFile=$2


python /home/kpradha1/projects/sanchari/filterCpg/filterCpg.py < $inFile > $outFile




