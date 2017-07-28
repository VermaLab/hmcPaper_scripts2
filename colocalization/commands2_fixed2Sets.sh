##########################################
#prepare the files for correlation with atac-seq


cat jd_mc_100.txt | cut -f1,2 | sed 's/.*/&\t+/' >| jd_mc_100.bed
cat jd_mc_lt50.txt | cut -f1,2 | sed 's/.*/&\t+/' >| jd_mc_lt50.bed
cat jd_mc_gt50.txt | cut -f1,2 | sed 's/.*/&\t+/' >| jd_mc_gt50.bed
cat jd_hmc_p05.txt | cut -f1,2 | sed 's/.*/&\t+/' >| jd_hmc_p05.bed

cat hpne_mc_100.txt | cut -f1,2 | sed 's/.*/&\t+/' >|  hpne_mc_100.bed
cat hpne_mc_lt50.txt | cut -f1,2 | sed 's/.*/&\t+/' >| hpne_mc_lt50.bed
cat hpne_mc_gt50.txt | cut -f1,2 | sed 's/.*/&\t+/' >| hpne_mc_gt50.bed
cat hpne_hmc_p05.txt | cut -f1,2 | sed 's/.*/&\t+/' >| hpne_hmc_p05.bed



qsub runPeakProfiles.sh jd_mc_100_jdAtac  jd_mc_100.bed /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000
qsub runPeakProfiles.sh jd_mc_lt50_jdAtac jd_mc_lt50.bed /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000
qsub runPeakProfiles.sh jd_mc_gt50_jdAtac jd_mc_gt50.bed /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000
qsub runPeakProfiles.sh jd_hmc_p05_jdAtac jd_hmc_p05.bed /home/kpradha1/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam 3000

qsub runPeakProfiles.sh hpne_mc_100_hpneAtac  hpne_mc_100.bed  /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000
qsub runPeakProfiles.sh hpne_mc_lt50_hpneAtac hpne_mc_lt50.bed /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000
qsub runPeakProfiles.sh hpne_mc_gt50_hpneAtac hpne_mc_gt50.bed /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000
qsub runPeakProfiles.sh hpne_hmc_p05_hpneAtac hpne_hmc_p05.bed /home/kpradha1/projects/sanchari/j455/hpne_sorted_mdup.bam 3000



python makeColocPlots.py jd_mc_100_jdAtac  
python makeColocPlots.py jd_mc_lt50_jdAtac 
python makeColocPlots.py jd_mc_gt50_jdAtac 
python makeColocPlots.py jd_hmc_p05_jdAtac 
python makeColocPlots.py hpne_mc_100_hpneAtac 
python makeColocPlots.py hpne_mc_lt50_hpneAtac
python makeColocPlots.py hpne_mc_gt50_hpneAtac
python makeColocPlots.py hpne_hmc_p05_hpneAtac



####################################
#now prepare the files for RNA correlation

#prepare the stranded position info
BAMJD=/home/kpradha1/projects/sanchari/rnaseq/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam
BAMHPNE=/home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam


qsub runPeakProfiles.sh jd_mc_100_jdRna  jd_mc_100.bed $BAMJD 3000
qsub runPeakProfiles.sh jd_mc_lt50_jdRna jd_mc_lt50.bed $BAMJD 3000
qsub runPeakProfiles.sh jd_mc_gt50_jdRna jd_mc_gt50.bed $BAMJD 3000
qsub runPeakProfiles.sh jd_hmc_p05_jdRna jd_hmc_p05.bed $BAMJD 3000

qsub runPeakProfiles.sh hpne_mc_100_hpneRna  hpne_mc_100.bed  $BAMHPNE 3000
qsub runPeakProfiles.sh hpne_mc_lt50_hpneRna hpne_mc_lt50.bed $BAMHPNE 3000
qsub runPeakProfiles.sh hpne_mc_gt50_hpneRna hpne_mc_gt50.bed $BAMHPNE 3000
qsub runPeakProfiles.sh hpne_hmc_p05_hpneRna hpne_hmc_p05.bed $BAMHPNE 3000



python makeColocPlots.py  jd_mc_100_jdRna  
python makeColocPlots.py jd_mc_lt50_jdRna 
python makeColocPlots.py jd_mc_gt50_jdRna 
python makeColocPlots.py jd_hmc_p05_jdRna 
python makeColocPlots.py  hpne_mc_100_hpneRna 
python makeColocPlots.py hpne_mc_lt50_hpneRna
python makeColocPlots.py hpne_mc_gt50_hpneRna
python makeColocPlots.py hpne_hmc_p05_hpneRna


##copy all the coloc.plots and rename them to a separate folder
mkdir pics2
mkdir pics3
ls hpne_{,h}mc_*/*.png jd_{,h}mc_*/*.png | sed 's/\(.*\)\/colocPlot\(.*\.png\)/cp & pics3\/\1\2/'






#look at JUST the profile and use full MC set, (not 400k max)

qsub runPeakProfiles_noHeatmap.sh jd_mc_lt50_jdRna jd_mc_lt50.bed $BAMJD 3000
qsub runPeakProfiles_noHeatmap.sh jd_mc_gt50_jdRna jd_mc_gt50.bed $BAMJD 3000

qsub runPeakProfiles_noHeatmap.sh hpne_mc_lt50_hpneRna hpne_mc_lt50.bed $BAMHPNE 3000
qsub runPeakProfiles_noHeatmap.sh hpne_mc_gt50_hpneRna hpne_mc_gt50.bed $BAMHPNE 3000



#########################################
#correlation with stat32
#have to strip off the "chr" from the positions 
cat jd_hmc_p05.bed | sed 's/^chr//' >| jd_hmc_p05_s3.bed
cat jd_mc_gt50.bed | sed 's/^chr//' >| jd_mc_gt50_s3.bed
cat hpne_hmc_p05.bed | sed 's/^chr//' >| hpne_hmc_p05_s3.bed
cat hpne_mc_gt50.bed | sed 's/^chr//' >| hpne_mc_gt50_s3.bed

qsub runPeakProfiles.sh jd_mc_gt50_jdStat3 jd_mc_gt50_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/jd13d-stat3.bam 3000
qsub runPeakProfiles.sh jd_hmc_p05_jdStat3 jd_hmc_p05_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/jd13d-stat3.bam 3000
qsub runPeakProfiles.sh hpne_mc_gt50_hpneStat3 hpne_mc_gt50_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/hpne-stat3.bam 3000
qsub runPeakProfiles.sh hpne_hmc_p05_hpneStat3 hpne_hmc_p05_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/hpne-stat3.bam 3000


qsub runPeakProfiles_noHeatmap.sh jd_mc_gt50_jdStat3 jd_mc_gt50_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/jd13d-stat3.bam 3000
qsub runPeakProfiles_noHeatmap.sh jd_hmc_p05_jdStat3 jd_hmc_p05_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/jd13d-stat3.bam 3000
qsub runPeakProfiles_noHeatmap.sh hpne_mc_gt50_hpneStat3 hpne_mc_gt50_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/hpne-stat3.bam 3000
qsub runPeakProfiles_noHeatmap.sh hpne_hmc_p05_hpneStat3 hpne_hmc_p05_s3.bed /home/kpradha1/projects/sanchari/stat3Bams/hpne-stat3.bam 3000



python makeColocPlots.py jd_mc_gt50_jdStat3  
python makeColocPlots.py jd_hmc_p05_jdStat3  
python makeColocPlots.py hpne_mc_gt50_hpneStat3  
python makeColocPlots.py hpne_hmc_p05_hpneStat3  





##############################3
####scratch







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
