

#atack seq peaks, RNA reads
########################3333
#run from cluster
#qsub runPeakProfiles.sh test2 hpneHmc_strandedPeakSummits.txt /home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam 3000
cut -f1,2,6 hpne_strandedPeakSummits.txt > hpne_atac.txt
cut -f1,2,6 jd_strandedPeakSummits.txt > jd_atac.txt


qsub runPeakProfiles.sh hpneAtacPeaks_hpneRna hpne_atac.txt /home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam 3000
qsub runPeakProfiles.sh jdAtacPeaks_jdRna jd_atac.txt /home/kpradha1/projects/sanchari/rnaseq/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam 3000


#need to prepare the atac-seq DIFF peaks
#use the macs2 run where the other samples was used as control

#1.  need to make the stranded peak summits files for the diffs
#.../sanchari/peaks_rna/annotationPeaks.r

cut -f1,2,6 hpne-jd_strandedPeakSummits.txt > hpne-jd_atac.txt
cut -f1,2,6 jd-hpne_strandedPeakSummits.txt > jd-hpne_atac.txt

#this keeps crashing.  maybe try to
#take out all the non chr* entries
grep -v hpne-jd_atac.txt -e "^chr."

grep -v jd-hpne_atac.txt -e "^chr[^_]\+_" > jd-hpne_atac1.txt
grep -v hpne-jd_atac.txt -e "^chr[^_]\+_" > hpne-jd_atac1.txt

qsub runPeakProfiles.sh hpne-jdAtacPeaks_hpneRna hpne-jd_atac1.txt /home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam 3000
qsub runPeakProfiles.sh jd-hpneAtacPeaks_jdRna jd-hpne_atac1.txt /home/kpradha1/projects/sanchari/rnaseq/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam 3000



#get plots of the JD gain/loss peak sets

#first, get stranded info for the peak sets
#/hpc_home/projects/sanchari/peaks_rna/annotatePeaks.r
#peaks.atac.jd.gain.1k.txt
#peaks.atac.jd.loss.1k.txt

cut -f1,2,6 peaks.atac.jd.gain.1k_stranded.txt > jdGain1k_atac.txt
cut -f1,2,6 peaks.atac.jd.loss.1k_stranded.txt > jdLoss1k_atac.txt

qsub runPeakProfiles.sh jdGain1kAtacPeaks_jdRna jdGain1k_atac.txt /home/kpradha1/projects/sanchari/rnaseq/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam 3000
qsub runPeakProfiles.sh jdLoss1kAtacPeaks_hpneRna jdLoss1k_atac.txt /home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam 3000





#hmc sites, RNA reads
#no..have to add strandedness info.
####################################

#todo
#.../sanchari/peaks_rna/annotationPeaks.r

qsub runPeakProfiles.sh hpneHmc_hpneRna hpne_strandedHmc.txt /home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam 3000
qsub runPeakProfiles.sh jdHmc_jdRna jd_strandedHmc.txt /home/kpradha1/projects/sanchari/rnaseq/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam 3000




#oxbs hmc sites, atac seq reads
###################################################33

#prep oxbs sites
cat hmc_sig0.01_commonSites1_jd.bed | cut -f1,2 | sed 's/.*/&\t+/' >| jd_hmc.txt
cat hmc_sig0.01_commonSites1_hpne.bed | cut -f1,2 | sed 's/.*/&\t+/' >| hpne_hmc.txt


qsub runPeakProfiles.sh hpneHmcSites_hpneAtac hpne_hmc.txt /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000
qsub runPeakProfiles.sh jdHmcSites_jdAtac jd_hmc.txt /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000




#oxbs mc sites, atac seq reads
###################################################33
#prep mc files
cat mc_commonSites1_jd.bed | cut -f1,2 | sed 's/.*/&\t+/' >| jd_mc.txt
cat mc_commonSites1_hpne.bed | cut -f1,2 | sed 's/.*/&\t+/' >| hpne_mc.txt

qsub runPeakProfiles.sh hpneMcSites_hpneAtac_400k hpne_mc.txt /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000
qsub runPeakProfiles.sh jdMcSites_jdAtac_400k jd_mc.txt /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000



#mc files with rna
################################
#run annotation to get strand info




#stat3 sites, atac seq reads
###########################################33

#prepare the stat3 peaks from bed files

cat jd13d_stat3_summits.bed | cut -f1,2 | sed 's/.*/&\t+/' | grep -v "hs37d5" | sed 's/.*/chr&/' >| jd_stat3.txt
cat hpne_stat3_summits.bed | cut -f1,2 | sed 's/.*/&\t+/'  | grep -v "hs37d5" | sed 's/.*/chr&/' >| hpne_stat3.txt


qsub runPeakProfiles.sh hpneStat3Sites_hpneAtac hpne_stat3.txt /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000
qsub runPeakProfiles.sh jdStat3Sites_jdAtac jd_stat3.txt /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000



#prepare stranded stat3 files with rna
#add chr to the front of each row
cat jd_stat3_summits_stranded.txt | sed 's/.*/chr&/' > jd_stat3_summits_stranded1.txt
cat hpne_stat3_summits_stranded.txt | sed 's/.*/chr&/' > hpne_stat3_summits_stranded1.txt


qsub runPeakProfiles.sh jdstat3Peaks_jdRna jd_stat3_summits_stranded1.txt /home/kpradha1/projects/sanchari/rnaseq/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam 3000
qsub runPeakProfiles.sh hpnestat3Peaks_hpneRna hpne_stat3_summits_stranded1.txt /home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam 3000



#move all the pics to a separate folder
ls */colocPlot.png | sed 's$\(.*\)\/\(.*\)$cp & pics1\/\1\.png$' | bash


########################################################################



bamfile1="/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"

samtools view -b $bamfile1 chr1:561524-567524 | samtools mpileup - 

samtools mpileup $bamfile1 -r chr1:561524-567524 | cut -f1-5


#run from my comp
python makeColocPlots.py test2 
