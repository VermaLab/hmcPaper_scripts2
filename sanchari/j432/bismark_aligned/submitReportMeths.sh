
#accidently name the bams .bam.bam
#fix that here
#ls *.bam.bam | sed 's/\(.*\.bam\)\.bam/mv & \1/'

#qsub reportMeth.sh J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted.bam
#qsub reportMeth.sh J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_merged_sorted.bam
#qsub reportMeth.sh J432_BHKVCGBCXX_Lane1_CTTGTA_JD-BS_merged_sorted.bam
#qsub reportMeth.sh J432_BHKVCGBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted.bam

#qsub sortBam.sh J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted_mdup_nsorted_clip.bam
#qsub sortBam.sh J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_merged_sorted_mdup_nsorted_clip.bam
#qsub sortBam.sh J432_BHKVCGBCXX_Lane1_CTTGTA_JD-BS_merged_sorted_mdup_nsorted_clip.bam
#qsub sortBam.sh J432_BHKVCGBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted_mdup_nsorted_clip.bam

qsub reportMeth.sh J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted_mdup_clip.bam
qsub reportMeth.sh J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_merged_sorted_mdup_clip.bam
qsub reportMeth.sh J432_BHKVCGBCXX_Lane1_CTTGTA_JD-BS_merged_sorted_mdup_clip.bam
qsub reportMeth.sh J432_BHKVCGBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted_mdup_clip.bam







