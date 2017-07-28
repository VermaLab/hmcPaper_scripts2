
ls ../j515/jd_sorted.bam
ls ../j455/jd_sorted.bam

#first remove duplicates
qsub markDuplicates.sh jd_sorted.bam


#then merge the bams
qsub mergeBams.sh jd_sorted_mdup_merged.bam ../j515/jd_sorted_mdup.bam ../j455/jd_sorted_mdup.bam


#now run peak finding 

qsub ~/projects/qsub_scripts/runMacs.sh jd_sorted_mdup_merged.bam

ln -s ../j455/hpne_sorted_mdup.bam .
qsub ~/projects/qsub_scripts/runMacs2.sh jd_sorted_mdup_merged.bam hpne_sorted_mdup.bam
qsub ~/projects/qsub_scripts/runMacs2.sh hpne_sorted_mdup.bam jd_sorted_mdup_merged.bam
