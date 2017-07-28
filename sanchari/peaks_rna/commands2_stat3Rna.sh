
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


