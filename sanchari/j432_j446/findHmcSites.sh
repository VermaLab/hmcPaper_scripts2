#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=1G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#   #$ -M kith.pradhan@einstein.yu.edu
module load python/2.7.8/gcc.4.4.7
module load scipy/0.14.0/python.2.7.8

oxFile=$1
oxbsFile=$2
outFile=$3


python findHmcSites.py $oxFile $oxbsFile $outFile




