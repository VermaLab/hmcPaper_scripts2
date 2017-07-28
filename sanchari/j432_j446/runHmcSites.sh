

./findHmcSites.r conversionCounts/HPNE_BS_22.txt.gz conversionCounts/HPNE_OxBS_22.txt.gz hmcSites_HPNE_22.txt




#for chr in {1..22} X Y
for chr in Y {22..1} X
do
    echo $chr
    ./findHmcSites.r conversionCounts/HPNE_BS_${chr}.txt.gz conversionCounts/HPNE_OxBS_${chr}.txt.gz hmcSites_HPNE_${chr}.txt
    ./findHmcSites.r conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz hmcSites_JD_${chr}.txt
done

for chr in Y {22..1} X
do
    echo $chr
    ./findHmcSites.r conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz hmcSites_JD_${chr}.txt
done

#run the stream based version of the prog
for chr in 1 2 3 
do
    echo $chr
    ./findHmcSites_stream.r conversionCounts/HPNE_BS_${chr}.txt.gz conversionCounts/HPNE_OxBS_${chr}.txt.gz stream_hmcSites_HPNE_${chr}.txt&
    ./findHmcSites_stream.r conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz stream_hmcSites_JD_${chr}.txt&
done


#run the python based version


for chr in Y {22..1} X
do
    echo $chr
    ./findHmcSites.py conversionCounts/HPNE_BS_${chr}.txt.gz conversionCounts/HPNE_OxBS_${chr}.txt.gz hmcSites_fixed/hmcSites_HPNE_${chr}.txt
    ./findHmcSites.py conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz hmcSites_fixed/hmcSites_JD_${chr}.txt
done


#12-12-16
#run on hpc cluster
#12-15-16
#reran, because a bug was found in findHmcSites.  Manmy sites were missing
#change script to only output sites with at least 4 reads in both bs and oxbs
for chr in Y {22..1} X
do
    echo $chr
    qsub findHmcSites.sh conversionCounts/HPNE_BS_${chr}.txt.gz conversionCounts/HPNE_OxBS_${chr}.txt.gz hmcSites_fixed2/hmcSites_HPNE_${chr}.txt
    qsub findHmcSites.sh conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz hmcSites_fixed2/hmcSites_JD_${chr}.txt
done

#why is chr17 so small in JD?
qsub findHmcSites.sh conversionCounts/JD_BS_17.txt.gz conversionCounts/JD_OxBs_17.txt.gz test1/hmcSites_JD_17.txt

#now filter them all for cpg sites
chr="Y"

qsub ../filterCpg/filterCpg.sh hmcSites_fixed/hmcSites_JD_${chr}.txt hmcSites_fixed/cpg/hmcSites_JD_${chr}_cpg.txt
qsub ../filterCpg/filterCpg.sh test1/hmcSites_JD_${chr}.txt test1/hmcSites_JD_${chr}_cpg.txt


#12-13-16
for chr in Y {22..1} X
do
    echo $chr
    qsub ../filterCpg/filterCpg.sh hmcSites_fixed2/hmcSites_JD_${chr}.txt hmcSites_fixed2/cpg/hmcSites_JD_${chr}_cpg.txt
    qsub ../filterCpg/filterCpg.sh hmcSites_fixed2/hmcSites_HPNE_${chr}.txt hmcSites_fixed2/cpg/hmcSites_HPNE_${chr}_cpg.txt
done


cd cpg



#join the jd and hpne sites into one
#run mergeSites.r


#combine all the chromsomes into one file
for chr in {1..22} X Y
do
    echo $chr
    cat hmcSites_BOTH_${chr}_cpg_common4.txt >> hmcSites_BOTH_all_cpg_common4.txt
done


