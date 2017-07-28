#!/bin/bash
#$ -cwd
#$ -j n
#$ -l h_vmem=10G
#$ -pe smp 1
#$ -S /bin/bash
#$ -m e
#usage:
#qsub runPeakProfiles.sh asdlk asdflka
#probably going to need more memory
module load scipy/0.14.0/python.2.7.8-atlas-3.11.30
module load python/2.7.8/gcc.4.4.7
module load HTSeq/0.6.1/python.2.7.8 
module load numpy/1.9.0/python.2.7.8-atlas-3.11.30
module load matplotlib/1.4.3/python.2.7.8



outputFolder="test3"
#siteFile="hpneHmc_strandedPeakSummits.txt"
siteFile="hpne_atac.txt"
#bamFile="/home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"
bamFile="/home/kpradhan/Desktop/hpc_home/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"
W="3000"

outputFolder=$1
siteFile=$2
bamFile=$3
W=$4


echo $outputFolder
echo $siteFile
echo $bamFile
echo $W


#make a folder for storing intermediate resulst
#otherwise, might run out of space on default /tmp
python getPeakProfiles_noHeatmap_preload.py $outputFolder $siteFile $bamFile $W

#add a $ after the # if want email notification at job completion
# -M kith.pradhan@einstein.yu.edu


