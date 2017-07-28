#step 1:  split each file into 40M line chunks
mkdir split
NL=40000000
for fq in `ls *.fq.gz`;
do
    nam=`basename $fq .gz`
    #get name without extension
    #echo $fq
    #echo $nam
    echo "split <(zcat $fq) -d -a 3 -l $NL split/${nam}.pt"
done >| split1.sh
chmod +x split1.sh


#step 2:  run bismark on the split files
#   and sort the resulting bam
mkdir bismark
cd bismark
for fq1 in `ls ../split/*1_val_1*`;
do
    fq2=`echo $fq1 | sed 's/\(.*\)1_val_1\(.*\)/\12_val_2\2/'`
    pname=`echo $fq1 | sed 's/.*\/\(.*\).1_val_1.fq.\(p.*\)/\1_\2/'`
    #echo $fq1
    #echo $fq2
    #echo $pname
    qsub ../runBismark.sh $pname $fq1 $fq2
done

#step 3.  merge the bams together
mkdir merged
cd merged
qsub ../mergeBams.sh J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted.bam ../bismark/*HPNE-1-BS_p*_sorted.bam
qsub ../mergeBams.sh J446_AHLT3HBCXX_Lane2_GCCAAT_HPNE-1-OxBS_merged_sorted.bam ../bismark/*HPNE-1-OxBS_p*_sorted.bam
qsub ../mergeBams.sh J446_BHLYGMBCXX_Lane1_CTTGTA_JD-BS_merged_sorted.bam ../bismark/*JD-BS_p*_sorted.bam
qsub ../mergeBams.sh J446_BHLYGMBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted.bam ../bismark/*JD-OxBS_p*_sorted.bam

#step 4.  mark duplicates
qsub ../markDuplicates.sh J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted.bam 
qsub ../markDuplicates.sh J446_AHLT3HBCXX_Lane2_GCCAAT_HPNE-1-OxBS_merged_sorted.bam
qsub ../markDuplicates.sh J446_BHLYGMBCXX_Lane1_CTTGTA_JD-BS_merged_sorted.bam
qsub ../markDuplicates.sh J446_BHLYGMBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted.bam 

#step 5.  clip overlap
qsub ../clipOverlaps.sh J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted_mdup.bam 
qsub ../clipOverlaps.sh J446_AHLT3HBCXX_Lane2_GCCAAT_HPNE-1-OxBS_merged_sorted_mdup.bam 
qsub ../clipOverlaps.sh J446_BHLYGMBCXX_Lane1_CTTGTA_JD-BS_merged_sorted_mdup.bam
qsub ../clipOverlaps.sh J446_BHLYGMBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted_mdup.bam 

#step 6.  methylation report generation
qsub ../reportMeth.sh J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted_mdup_clip.bam    
qsub ../reportMeth.sh J446_AHLT3HBCXX_Lane2_GCCAAT_HPNE-1-OxBS_merged_sorted_mdup_clip.bam  
qsub ../reportMeth.sh J446_BHLYGMBCXX_Lane1_CTTGTA_JD-BS_merged_sorted_mdup_clip.bam        
qsub ../reportMeth.sh J446_BHLYGMBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted_mdup_clip.bam      



