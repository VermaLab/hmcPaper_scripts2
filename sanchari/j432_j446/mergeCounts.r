list.files(pattern="*.bedGraph")

list.files(pattern="*HPNE-1[-_]BS*")

chrom = "chr10"

x1 = read.table(pipe("grep J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted_mdup_clip.bedGraph -e \"^chr17\\s\" | cut -f1,2,5,6"))










