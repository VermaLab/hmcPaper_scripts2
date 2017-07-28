ngsplotdb.py list

tss, tes, genebody, exon, cgi, enhancer, dhs
testbam=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam
ngs.plot.r -G hg19 -R tss -C $testbam -O test1plot_tss -T testplot -L 3000 
ngs.plot.r -G hg19 -R genebody -C $testbam -O test1plot_genebody -T testplot -L 3000 
ngs.plot.r -G hg19 -R exon -C $testbam -O test1plot_exon -T testplot -L 3000 
ngs.plot.r -G hg19 -R enhancer -C $testbam -O test1plot_enhancer -T testplot -L 3000 
ngs.plot.r -G hg19 -R dhs -C $testbam -O test1plot_dhs -T testplot -L 3000 

#make a test bed file
#chr start stop
ls ../../

#prepare the hmc site bed files for jd and hpne
tail -n +2 ../j432_j446/hmcSites/hmcSites_HPNE_all.txt | sed 's/"//g' | awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| sites_hpne_all.bed
tail -n +2 ../j432_j446/hmcSites/hmcSites_JD_all.txt | sed 's/"//g' | awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| sites_jd_all.bed

#prepare a list of random sites



#some smaller single chrom beds
tail -n +2 ../j432_j446/hmcSites/hmcSites_HPNE_10.txt | sed 's/"//g' | awk '{print $1,"\t",$2-2000,"\t",$2+2000}' | sed 's/ //'g >| s10.bed
tail -n +2 ../j432_j446/hmcSites/hmcSites_HPNE_10.txt | sed 's/"//g' | awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| s10.bed
time ngs.plot.r -G hg19 -R bed -E s10.bed -C $bam -O test2plot_bed -T testplot -P 6 -L 3000 

#prepare bed files for chr19
tail -n +2 ../j432_j446/hmcSites/hmcSites_HPNE_all.txt | sed 's/"//g' | awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g | grep -e "chr19" >| sites_hpne_19.bed
tail -n +2 ../j432_j446/hmcSites/hmcSites_JD_all.txt | sed 's/"//g' | awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g | grep -e "chr19" >| sites_jd_19.bed


bam1=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/hpne_sorted_mdup.bam
#bam2=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/jd_sorted_mdup.bam
bam2=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam
#time ngs.plot.r -G hg19 -R bed -E sites_hpne_19.bed -C $bam1 -O hpne_19 -T hpne_19 -P 6 -L 3000 
#cat test1_jdbad.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' > test2_jdbad.txt
#cat test1_jdmc.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' > test2_jdmc.txt
#cat test1_jdmc_top10k.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' > test2_jdmc_top10k.txt
#cat test1_jdmc_top20k-10k.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' >| test2_jdmc_top20k-10k.txt
#cat test1_jdmc_top30k-10k.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' >| test2_jdmc_top30k-10k.txt
cat test1_jdmc_top60k-10k.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' >| test2_jdmc_top60k-10k.txt
cat test1_jdmc_top100k.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' >| test2_jdmc_top100k.txt
cut -f1-3 test2_jdmc_top30k-10k.txt > test3_jdmc_top30k-10k.txt
sed test2_jdmc_top30k-10k.txt -e 's/+/-/' > test4_jdmc_top30k-10k.txt
time ngs.plot.r -G hg19 -R bed -E sites_jd_19.bed -C $bam2 -O jd_19 -T jd_19 -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdhmc.txt -C $bam2 -O jd_test1 -T jd_test -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdmc.txt -C $bam2 -O jd_mc_test1 -T jd_mc_test -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdmc_top10k.txt -C $bam2 -O jdmc_top10k.txt -T jdmc_top10k.txt -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdmc_top20k-10k.txt -C $bam2 -O top20k-10k -T jdmc_top20k-10k.txt -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test3_jdmc_top30k-10k.txt -C $bam2 -O top30k-10k_v2 -T jdmc_top30k-10k_v2 -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test4_jdmc_top30k-10k.txt -C $bam2 -O top30k-10k_v4 -T jdmc_top30k-10k_v4 -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdmc_top60k-10k.txt -C $bam2 -O jdmc_top60k-10k -T jdmc_top60k-10k -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdmc_top100k.txt -C $bam2 -O jdmc_100k -T jdmc_100k -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -E test2_jdbad.txt -C $bam2 -O jdbad_test1 -T jdbad_test -P 6 -L 3000 

time ngs.plot.r -G hg19 -R bed -E sites_rand_chr19.bed -C config1.txt -O jd_rand19 -T jd_rand19 -P 6 -L 3000 
time ngs.plot.r -G hg19 -R bed -C config2.txt -O jd_batch -T jd_batch -P 6 -L 3000 

#prepare the hmcSite bed file from the common4 hmcSites.
#TODO:  messed up JD common sites bed..must redo
cat ../j432_j446/hmcSites/common4_hmcSites_JD.txt | cut -f1,2 | sed 's/#//' | \
    awk '{print $1,"\t",$2-1000,"\t",$2+1000}' | sed 's/ //'g >| commonSites_jd.bed
cat ../j432_j446/hmcSites/common4_hmcSites_HPNE.txt | cut -f1,2 | sed 's/#//' | \
    awk '{print $1,"\t",$2-1000,"\t",$2+1000}' | sed 's/ //'g >| commonSites_hpne.bed


#dont' forget to take out spaces.
#ngsplot only works with tabs
cat ../j432_j446/hmcSites/common4_hmcSites_JD.txt | cut -f1,2 | sed 's/#//' | \
    awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| commonSites1_jd.bed
cat ../j432_j446/hmcSites/common4_hmcSites_HPNE.txt | cut -f1,2 | sed 's/#//' | \
    awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| commonSites1_hpne.bed



time ngs.plot.r -S 0.1 -G hg19 -R bed -E commonSites1_hpne.bed -C $bam1 -O hpne_all -T hpne_all -P 2 -L 3000
#time ngs.plot.r -S 0.01 -G hg19 -R bed -E commonSites1_hpne.bed -C $bam2 -O jdBam_hpneSites_all -T jdBam_hpneSites_0.01_all -P 6 -L 3000
time ngs.plot.r -S 0.1 -G hg19 -R bed -E commonSites1_jd.bed -C $bam2 -O jd_all -T jd_all -P 2 -L 3000 

time ngs.plot.r -S 0.01 -G hg19 -R bed -C config1.txt -O atacSeqBatch -T atacSeqBatch -P 6 -L 3000 



bam1=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/hpne_sorted_mdup.bam
bam2=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam
#try it with a more significant set of hmc values
cat ../j432_j446/hmcSites/common4_hmcSites_JD.txt | cut -f1,2,8 | sed 's/#//' | \
    awk '{if ($3 < 0.01) print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| hmc_sig0.01_commonSites1_jd.bed
cat ../j432_j446/hmcSites/common4_hmcSites_HPNE.txt | cut -f1,2,8 | sed 's/#//' | \
    awk '{if ($3 < 0.01) print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| hmc_sig0.01_commonSites1_hpne.bed

time ngs.plot.r  -G hg19 -R bed -E hmc_sig0.01_commonSites1_hpne.bed -C $bam1 -O hmc_sig0.01_hpne_all -T hmc_sig0.01_hpne_all -P 2 -L 3000
time ngs.plot.r  -G hg19 -R bed -E hmc_sig0.01_commonSites1_jd.bed -C $bam2 -O hmc_sig0.01_jd_all -T hmc_sig0.01_jd_all -P 2 -L 3000 


#prepare bed files for the stat3 peaks
bam1=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/hpne_sorted_mdup.bam
bam2=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam

time ngs.plot.r  -G hg19 -R bed -E hpne_stat3_summits.bed -C $bam1 -O hpne_stat3Peaks_atac -T hpne_stat3Peaks_atac -P 2 -L 3000
time ngs.plot.r  -G hg19 -R bed -E jd13d_stat3_summits.bed -C $bam2 -O jd_stat3Peaks_atac -T jd_stat3Peaks_atac -P 2 -L 3000



#make plots for the MC sites
###############################################


#dont' forget to take out spaces.
#ngsplot only works with tabs
cat ../j432_j446/hpc_mcSites/common4_mcSites_JD.txt | cut -f1,2 | sed 's/#//' | \
    awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| mc_commonSites1_jd.bed
cat ../j432_j446/hpc_mcSites/common4_mcSites_HPNE.txt | cut -f1,2 | sed 's/#//' | \
    awk '{print $1,"\t",$2,"\t",$2+1}' | sed 's/ //'g >| mc_commonSites1_hpne.bed


bam1=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/hpne_sorted_mdup.bam
bam2=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam

time ngs.plot.r -S 0.05 -G hg19 -R bed -E mc_commonSites1_hpne.bed -C $bam1 -O mc_hpne_0.05_all -T mc_hpne_0.05_all -P 2 -L 3000
time ngs.plot.r -S 0.05 -G hg19 -R bed -E mc_commonSites1_jd.bed -C $bam2 -O mc_jd_0.05_all -T mc_jd_0.05_all -P 2 -L 3000 

time ngs.plot.r -S 0.01 -G hg19 -R bed -E mc_commonSites1_jd.bed -C $bam2 -O mc_jd_0.01_all -T mc_jd_0.01_all -P 3 -L 3000 
time ngs.plot.r -S 0.01 -G hg19 -R bed -E mc_commonSites1_hpne.bed -C $bam1 -O mc_hpne_0.01_all -T mc_hpne_0.01_all -P 3 -L 3000

time ngs.plot.r -S 0.1 -G hg19 -R bed -E commonSites1_jd.bed -C $bam2 -O hmc_jd_0.1_all -T hmc_jd_0.1_all -P 3 -L 3000 
time ngs.plot.r -S 0.1 -G hg19 -R bed -E commonSites1_hpne.bed -C $bam1 -O hmc_hpne_0.1_all -T hmc_hpne_0.1_all -P 3 -L 3000




#try to do the same with RNA seq data

bam1=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam
bam2=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam

samtools index $bam1
samtools index $bam2

time ngs.plot.r -S 0.01 -G hg19 -R bed -E mc_commonSites1_hpne.bed -C $bam1 -O mc_rna_hpne_0.01_all -T mc_rna_hpne_0.01_all -P 2 -L 3000
time ngs.plot.r -S 0.01 -G hg19 -R bed -E mc_commonSites1_jd.bed -C $bam2 -O mc_rna_jd_0.01_all -T mc_rna_jd_0.01_all -P 2 -L 3000 
time ngs.plot.r -RR 300 -S 0.05 -G hg19 -R bed -E mc_commonSites1_hpne.bed -C $bam1 -O mc_rna_hpne_0.05_all -T mc_rna_hpne_0.05_all -P 2 -L 3000
time ngs.plot.r -RR 300 -S 0.05 -G hg19 -R bed -E mc_commonSites1_jd.bed -C $bam2 -O mc_rna_jd_0.05_all -T mc_rna_jd_0.05_all -P 2 -L 3000 





#filter out mitichondria 
samtools idxstats input.bam | cut -f 1 | grep -v MT | xargs samtools view -b input.bam > output.bam

#get the numbe of good reads

#3844
#read unmapped
#not primary alignment
#read fails platform/vendor quality checks
#read is PCR or optical duplicate
#supplementary alignment
#64
#first in pair
samtools view -f 64 -F 3844 -c $bam1
samtools view -f 64 -F 3844 -c $bam1 chrM
samtools view -f 64 -F 3844 -c $bam2
samtools view -f 64 -F 3844 -c $bam2 chrM






