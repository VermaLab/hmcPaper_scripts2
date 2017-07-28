#get a list of differential peaks between jd and hpne.
#the peaks in one but not the other.

#see if stat3 appears in either of the differential peaks

#want to find transciption factor bindings at these diff peaks
#MOTIF analysis?

#rna coverage plot at these diff peak summits

#prepare peak strand info for the hmc sites





#RNA bams
bam1=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam
bam2=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam

#atacseq peaks
peak1=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/macs_hpne_sorted_mdup/hpne_sorted_mdup_summits.bed
peak2=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_jd_sorted_mdup_merged/jd_sorted_mdup_merged_summits.bed

#atacseq peaks with strand information
peak1=hpne_strandedPeakSummits.txt
peak2=jd_strandedPeakSummits.txt

#hmc peaks with strand
peak1=hpneHmc_strandedPeakSummits_0.001.txt
peak2=jdHmc_strandedPeakSummits_0.001.txt

cat $peak1 | awk -v OFS='\t' '{print $1,$2,$3,"1234","5678",$4}' >| p1.txt
cat $peak2 | awk -v OFS='\t' '{print $1,$2,$3,"1234","5678",$4}' >| p2.txt
hpneHmc_strandedPeakSummits.txt


cat commonSites1_hpne.sorted.bed | awk -v OFS='\t' '{print $1,$2,$3,"1234","5678",$4}' >| p3.txt

peak1=hpneHmc_strandedPeakSummits_0.001.txt
peak2=jdHmc_strandedPeakSummits_0.001.txt

peak1=commonSites1_hpne.sorted.bed
peak2=commonSites1_jd.sorted.bed

time ngs.plot.r -SS same -G hg19 -R bed -E p1.txt -C $bam1 -O test1 -T test1 -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E p2.txt -C $bam2 -O test2 -T test2 -P 2 -L 3000

time ngs.plot.r -S 0.1 -SS same -G hg19 -R bed -E p3_stranded.txt -C $bam1 -O test3 -T test3 -P 2 -L 3000

time ngs.plot.r -SS same -G hg19 -R bed -E $peak1 -C $bam1 -O hpne_strandedHmcPeaks_fisher0.001_strandedRna -T hpne_strandedHmcPeaks_fisher0.001_strandedRna -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E $peak2 -C $bam2 -O jd_strandedHmcPeaks_fisher0.001_strandedRna -T jd_strandedHmcPeaks_fisher0.001_strandedRna -P 2 -L 3000


time ngs.plot.r -S 0.1 -SS same -G hg19 -R bed -E $peak1 -C $bam1 -O hpne_strandedHmcPeaks_0.01_strandedRna -T hpne_strandedHmcPeaks_0.1_strandedRna -P 2 -L 3000
time ngs.plot.r -S 0.1 -SS same -G hg19 -R bed -E $peak2 -C $bam2 -O jd_strandedHmcPeaks_0.01_strandedRna -T jd_strandedHmcPeaks_0.1_strandedRna -P 2 -L 3000


#hmc use same stranded reads for the stranded peaks
time ngs.plot.r -SS same -G hg19 -R bed -E $peak1 -C $bam1 -O hpneHmc_strandedAtacPeaks_strandedRna -T hpneHmc_strandedAtacPeaks_strandedRna -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E $peak2 -C $bam2 -O jdHmc_strandedAtacPeaks_strandedRna -T jdHmc_strandedAtacPeaks_strandedRna -P 2 -L 3000


#use same stranded reads for the stranded peaks
time ngs.plot.r -SS same -G hg19 -R bed -E $peak1 -C $bam1 -O hpne_strandedAtacPeaks_strandedRna -T hpne_strandedAtacPeaks_strandedRna -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E $peak2 -C $bam2 -O jd_strandedAtacPeaks_strandedRna -T jd_strandedAtacPeaks_strandedRna -P 2 -L 3000



#time ngs.plot.r -G hg19 -R bed -E $peak1 -C $bam1 -O atac_rna_hpne_all -T atac_rna_hpne_all -P 2 -L 3000
#time ngs.plot.r -G hg19 -R bed -E $peak2 -C $bam2 -O atac_rna_jd_all -T atac_rna_jd_all -P 2 -L 3000 

time ngs.plot.r -G hg19 -R bed -E peaks_test1.bed -C $bam1 -O testatac_rna_hpne_all -T testatac_rna_hpne_all -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E peaks_test1.bed -C $bam1 -O test2atac_rna_hpne_all -T test2atac_rna_hpne_all -P 2 -L 3000


#####################
#new stuff
#



#try it with hmc peaks?
time ngs.plot.r -G hg19 -R bed -E jdHmc_strandedPeakSummits.txt -C $bam2 -O jdHmc0.01_jdRna -T jdHmc0.01_jdRna -P 2 -L 3000
time ngs.plot.r -G hg19 -R bed -E hpneHmc_strandedPeakSummits.txt -C $bam1 -O hpneHmc0.01_hpneRna -T hpneHmc0.01_hpneRna -P 2 -L 3000
sort -k1,1 -k2,2n jdHmc_strandedPeakSummits.txt | ../filterCpg/filterCpg.py | cut -f1-4 >| jdHmc_cpg_strandedPeakSummits.txt
sort -k1,1 -k2,2n hpneHmc_strandedPeakSummits.txt | ../filterCpg/filterCpg.py | cut -f1-4 >| hpneHmc_cpg_strandedPeakSummits.txt
time ngs.plot.r -SS same -G hg19 -R bed -E jdHmc_cpg_strandedPeakSummits.txt -C $bam2 -O jdHmc0.01cpg_jdRna -T jdHmc0.01cpg_jdRna -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E hpneHmc_cpg_strandedPeakSummits.txt -C $bam1 -O hpneHmc0.01cpg_hpneRna -T hpneHmc0.01cpg_hpneRna -P 2 -L 3000

#TODO:
#instead of the 0.01 hmc sites...try will all of the hmc sites
#but use a random shuf of them.  


#try atacseq diff peaks and RNA
time ngs.plot.r -G hg19 -R bed -E jd-hpne_strandedPeakSummits.txt -C $bam2 -O jd-hpneAtac_jdRna -T jd-hpneAtac_jdRna -P 2 -L 3000
time ngs.plot.r -G hg19 -R bed -E hpne-jd_strandedPeakSummits.txt -C $bam1 -O hpne-jdAtac_jdRna -T hpne-jdAtac_jdRna -P 2 -L 3000


#TODO
#try it with mc peaks?


#scratch

hpneHmc_strandedPeakSummits.txt
cat jdHmc_strandedPeakSummits.txt | awk -v OFS='\t' '{print $1,$2,$3,"1234","5678",$4}' >| p1.txt
cat $peak2 | awk -v OFS='\t' '{print $1,$2,$3,"1234","5678",$4}' >| p2.txt

