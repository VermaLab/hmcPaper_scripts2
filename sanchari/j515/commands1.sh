
qsub runBowtie.sh jd  J515_BHMJHHBCXX_Lane2_TCCTGA_JD13D.1_val_1.fq.gz J515_BHMJHHBCXX_Lane2_TCCTGA_JD13D.2_val_2.fq.gz


#first remove duplicates
qsub markDuplicates.sh jd_sorted.bam



#run peak finding on this sample alone



qsub ~/projects/qsub_scripts/runMacs.sh jd_sorted_mdup.bam
#can't have .. in file names.
#make a link instead
#qsub ~/projects/qsub_scripts/runMacs2.sh jd_sorted_mdup.bam ../j455/hpne_sorted_mdup.bam

ln -s ../j455/hpne_sorted_mdup.bam .
qsub ~/projects/qsub_scripts/runMacs2.sh jd_sorted_mdup.bam hpne_sorted_mdup.bam


