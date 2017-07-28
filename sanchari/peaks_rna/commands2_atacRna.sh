
#RNA bams
bam1=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam
bam2=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam



#atacseq peaks with strand information
peak1=peaks.atac.jd.gain.1k_stranded.txt
peak2=peaks.atac.jd.loss.1k_stranded.txt



time ngs.plot.r -SS same -G hg19 -R bed -E $peak1 -C $bam1 -O peaks.atac.jd.gain.1k -T peaks.atac.jd.gain.1k -P 2 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E $peak2 -C $bam2 -O peaks.atac.jd.loss.1k -T peaks.atac.jd.loss.1k -P 2 -L 3000

python makeColocPlots.py jdLoss1kAtacPeaks_hpneRna
python makeColocPlots.py jdGain1kAtacPeaks_jdRna


#test out the new atac peaks
bam1=/home/kpradhan/Desktop/data/sanchari/pancreatic/rnaseq/bam/JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam
cut -f1,2,4 test1_jdmc_top100k_stranded.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' > test2_jdmc_top100k_stranded.txt
cut -f1,2,4 test1_jdhmc_stranded.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' > test2_jdhmc_stranded.txt

#cat test1_jdmc.txt | awk '{printf("%s\t%d\t%d\t%s\n", $1, $2, $2+1, $3)}' > test2_jdmc.txt
time ngs.plot.r -SS same -G hg19 -R bed -E test2_jdmc_top100k_stranded.txt -C $bam1 -O test_jdmc_top100k -T test_jdmc_top100k -P 6 -L 3000
time ngs.plot.r -SS same -G hg19 -R bed -E test2_jdhmc_stranded.txt -C $bam1 -O test_jdhmc -T test_jdhmc -P 6 -L 3000

