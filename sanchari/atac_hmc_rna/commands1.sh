
bam1=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam
bam2=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam



#atacseq peaks with strand information
#peak1=jd_strandedAtacPeakSummits.txt
#peak2=hpne_strandedAtacPeakSummits.txt


#all
#time ngs.plot.r -SS same -G hg19 -R bed -E $peak1 -C $bam1 -O peaks.atac.jd -T peaks.atac.jd -P 2 -L 3000
#time ngs.plot.r -SS same -G hg19 -R bed -E $peak2 -C $bam2 -O peaks.atac.hpne -T peaks.atac.hpne -P 2 -L 3000

#separated by hmc
time ngs.plot.r -SS same -G hg19 -R bed -E jd_strandedAtacPeaks_hmc.txt -C $bam1 -O jd_strandedAtacPeaks_hmc -T peaks.atac.jd.hmc -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E jd_strandedAtacPeaks_nohmc.txt -C $bam1 -O jd_strandedAtacPeaks_nohmc -T peaks.atac.jd.nohmc -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E hpne_strandedAtacPeaks_hmc.txt -C $bam2 -O hpne_strandedAtacPeaks_hmc -T peaks.atac.hpne.hmc -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E hpne_strandedAtacPeaks_nohmc.txt -C $bam2 -O hpne_strandedAtacPeaks_nohmc -T peaks.atac.hpne.nohmc -P 2 -L 3000
