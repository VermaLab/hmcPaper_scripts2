#testing
tail -50 J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted_mdup_clip.bedGraph >> file1.txt
tail -50 J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted_mdup_clip.bedGraph >> file2.txt
./mergeCounts.pl <(cut -f1,2,5,6 file1.txt) <(cut -f1,2,5,6 file2.txt) chr10



function mergeCounts {
    file1=$1
    file2=$2
    outfile=$3

    for chr in {1..22} X Y
    do
        echo $chr
        time ./mergeCounts.pl <(cut -f1,2,5,6 $file1 | grep -e "^chr$chr\s") <(cut -f1,2,5,6 $file2 | grep -e "^chr$chr\s") chr$chr | gzip > ${outfile}_${chr}.txt.gz
    done
}

mergeCounts\
    file1.txt\
    file2.txt\
    test

mergeCounts\
    J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted_mdup_clip.bedGraph\
    J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted_mdup_clip.bedGraph\
    HPNE_BS

mergeCounts\
    J432_AHHH5NBCXX_Lane2_GCCAAT_HPNE-1OxBS_merged_sorted_mdup_clip.bedGraph\
    J446_AHLT3HBCXX_Lane2_GCCAAT_HPNE-1-OxBS_merged_sorted_mdup_clip.bedGraph\
    HPNE_OxBS

mergeCounts\
    J432_BHKVCGBCXX_Lane1_CTTGTA_JD-BS_merged_sorted_mdup_clip.bedGraph\
    J446_BHLYGMBCXX_Lane1_CTTGTA_JD-BS_merged_sorted_mdup_clip.bedGraph\
    JD_BS

mergeCounts\
    J432_BHKVCGBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted_mdup_clip.bedGraph\
    J446_BHLYGMBCXX_Lane2_GTGAAA_JD-OxBS_merged_sorted_mdup_clip.bedGraph\
    JD_OxBs





#########################################

file1=J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS_merged_sorted_mdup_clip.bedGraph
file2=J446_AHLT3HBCXX_Lane1_ACAGTG_HPNE-1-BS_merged_sorted_mdup_clip.bedGraph
outfile=HPNE_BS
chr=22
echo $chr
time ./mergeCounts.pl <(cut -f1,2,5,6 $file1 | grep -e "^chr$chr\s") <(cut -f1,2,5,6 $file2 | grep -e "^chr$chr\s") chr$chr | gzip > ${outfile}_${chr}.txt.gz



#file1=file1.txt
#file2=file2.txt
#outfile=test
for chr in {1..22} X Y
do
    echo $chr
    time ./mergeCounts.pl <(cut -f1,2,5,6 $file1 | grep -e "^chr$chr\s") <(cut -f1,2,5,6 $file2 | grep -e "^chr$chr\s") chr$chr | gzip > ${outfile}_${chr}.txt.gz
done



#create set of "common sites".  
#sites that have at least 4 reads in all samples
ls conversionCounts/JD_BS_*.txt.gz
zcat conversionCounts/JD_BS_*.txt.gz | ./filterSitesByCount.pl 4 >| commonSites/common4_JD_BS.txt
zcat conversionCounts/HPNE_BS_*.txt.gz | ./filterSitesByCount.pl 4 >| commonSites/common4_HPNE_BS.txt
zcat conversionCounts/JD_OxBs_*.txt.gz | ./filterSitesByCount.pl 4 >| commonSites/common4_JD_OxBS.txt
zcat conversionCounts/HPNE_OxBS_*.txt.gz | ./filterSitesByCount.pl 4 > commonSites/common4_HPNE_OxBS.txt

#use my own computer to do the sorting
time sort --parallel=6 -T tmp commonSites/common4_JD_BS.txt >| commonSites/sorted_common4_JD_BS.txt
time sort --parallel=6 -T tmp commonSites/common4_HPNE_BS.txt >| commonSites/sorted_common4_HPNE_BS.txt
time sort --parallel=6 -T tmp commonSites/common4_JD_OxBS.txt >| commonSites/sorted_common4_JD_OxBS.txt
time sort --parallel=6 -T tmp commonSites/common4_HPNE_OxBS.txt >| commonSites/sorted_common4_HPNE_OxBS.txt


#join all the files, retaining only those present in all
cd commonSites
time join -t# sorted_common4_JD_OxBS.txt <(join -t# sorted_common4_JD_BS.txt <(join -t# sorted_common4_HPNE_BS.txt sorted_common4_HPNE_OxBS.txt)) >| common4.txt

#get the common sites only looking at a single sample
cd commonSites
time join -t# sorted_common4_JD_OxBS.txt sorted_common4_JD_BS.txt  >| common4_JD.txt&
time join -t# sorted_common4_HPNE_OxBS.txt sorted_common4_HPNE_BS.txt  >| common4_HPNE.txt&


#$ wc -l common4.txt
#527710019 common4.txt


#now lets filter away the hmcsites that aren't common.
#1.  sort hmcsites by chr pos alphabetically
#2.  join with commons sites
time join -t# ../commonSites/common4.txt <(tail -n +2 hmcSites_HPNE_all.txt | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/') >| common4_hmcSites_HPNE.txt
time join -t# ../commonSites/common4.txt <(tail -n +2 hmcSites_JD_all.txt | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/') >| common4_hmcSites_JD.txt

#wc -l common4_hmcSites_*
#  4657389 common4_hmcSites_HPNE.txt
#  6067629 common4_hmcSites_JD.txt
# 10725018 total



###############
#scratch
time join -t# ../commonSites/common4.txt <(head ../commonSites/common4.txt )
time join -t# ../commonSites/common4.txt <(tail -n +2 hmcSites_HPNE_all.txt | head -5000 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' | sort -k 1 -t#) 
time join -t# ../commonSites/common4.txt <(tail -n +2 hmcSites_HPNE_all.txt | grep -e "chr10" | head -5000 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' | sort -k 1 -t#) 

tail -n +2 hmcSites_HPNE_all.txt | grep -e "chr10" | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' | sort -k 1 -t# >| test1.txt
tail -n +2 hmcSites_HPNE_all.txt | grep -e "chr10" | sort -k 1,2  >| test1a.txt
tail -n +2 hmcSites_HPNE_all.txt | grep -e "chr10" | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' >| test1b.txt

join -t# ../commonSites/common4.txt test1.txt > test2.txt
#sort by chr then pos numerically



