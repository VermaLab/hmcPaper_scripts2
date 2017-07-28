qsub runBowtie.sh hpne J455_BC8P19ACXX_Lane4_AGGCAGAA-CTCTCTAT_HPNE.1_val_1.fq.gz J455_BC8P19ACXX_Lane4_AGGCAGAA-CTCTCTAT_HPNE.2_val_2.fq.gz

qsub runBowtie.sh jd  J455_BC8P19ACXX_Lane4_TCCTGAGC-AGAGTAGA_JD13D.1_val_1.fq.gz J455_BC8P19ACXX_Lane4_TCCTGAGC-AGAGTAGA_JD13D.2_val_2.fq.gz



SNAME=hpne
FQ1=J455_BC8P19ACXX_Lane4_AGGCAGAA-CTCTCTAT_HPNE.1_val_1.fq.gz
FQ2=J455_BC8P19ACXX_Lane4_AGGCAGAA-CTCTCTAT_HPNE.2_val_2.fq.gz
REF=../../hg19/hg19

qsub runMacs.sh hpne_sorted_mdup.bam
qsub runMacs.sh jd_sorted_mdup.bam



qsub runMacs2.sh jd_sorted_mdup.bam hpne_sorted_mdup.bam
qsub runMacs2.sh hpne_sorted_mdup.bam jd_sorted_mdup.bam


