#prepare files for mc discovery.

#both BS and OxBS should have no converted bases
#1.  write an r script that goes through each chromosome
#   writes the site position if both Bs and oxbs have 0% conversino

#run the stream based version of the prog
for chr in Y {22..1} X
do
    #echo $chr
    echo "./findMcSites_stream.r conversionCounts/HPNE_BS_${chr}.txt.gz conversionCounts/HPNE_OxBS_${chr}.txt.gz mcSites/stream_mcSites_HPNE_${chr}.txt"
    echo "./findMcSites_stream.r conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz mcSites/stream_mcSites_JD_${chr}.txt"
done | parallel -j 6


#messed up conversinoCountsHPNE_BS_22
#either redo....or skip it and take it from results above.


#do the same using the hpc
chr=21
echo $chr
qsub runFindMcSites.sh test1/HPNE_BS_${chr}.txt.gz test1/HPNE_OxBS_${chr}.txt.gz hpc_mcSites/stream_mcSites_HPNE_${chr}.txt

#editted this 12-12-16
for chr in Y {22..1} X
do
    echo $chr
    qsub runFindMcSites.sh conversionCounts/HPNE_BS_${chr}.txt.gz conversionCounts/HPNE_OxBS_${chr}.txt.gz mcSites_fixed/mcSites_HPNE_${chr}.txt
    qsub runFindMcSites.sh conversionCounts/JD_BS_${chr}.txt.gz conversionCounts/JD_OxBs_${chr}.txt.gz mcSites/mcSites_JD_${chr}.txt
done


#ran this 12-12-16
#the new faster way using unix tools
for chr in Y {22..1} X
do
    echo $chr
    time join <(zcat conversionCounts/HPNE_BS_${chr}.txt.gz | cut -f2,3,4 | awk '{ printf "%010i %i %i\n" , $1 , $2 , $3 }') <(zcat conversionCounts/HPNE_OxBS_${chr}.txt.gz | cut -f2,3,4 | awk '{ printf "%010i %i %i\n" , $1 , $2 , $3 }') | awk '{if ($2 == $3 && $3 >= 4 && $4 == $5 && $5 >= 4) print $0}' | sed 's/^0*//' | sed -e "s/.*/chr$chr &/" > mcSites_fixed/mcSites_HPNE_${chr}.txt
    time join <(zcat conversionCounts/JD_BS_${chr}.txt.gz | cut -f2,3,4 | awk '{ printf "%010i %i %i\n" , $1 , $2 , $3 }') <(zcat conversionCounts/JD_OxBs_${chr}.txt.gz | cut -f2,3,4 | awk '{ printf "%010i %i %i\n" , $1 , $2 , $3 }') | awk '{if ($2 == $3 && $3 >= 4 && $4 == $5 && $5 >= 4) print $0}' | sed 's/^0*//' | sed -e "s/.*/chr$chr &/" > mcSites_fixed/mcSites_JD_${chr}.txt
done

#join the BS and OXBS files
    #unzip
    #take out the chr* field
    #pad the position with zeros so it's lexigraphical sorted
#retain the rows that have 0 conversion and at least 4 reads in both
#unpad the position
#tack on the chr* field

#merge all the chromos into one file per sample
cd hpc_mcSites
for chr in {1..22} X Y
do
    echo $chr
    cat stream_mcSites_HPNE_${chr}.txt >> mcSites_HPNE_all.txt
    cat stream_mcSites_JD_${chr}.txt >> mcSites_JD_all.txt
done





#join with the list of common4 sites
#1.  sort mcsites by chr pos alphabetically
#       do this by adding a # between the 2nd and 3rd columns
#       effectively merging chr/pos into 1 field, and rest into another field
#    turn spaces into tabs
#2.  join with commons sites
time sort -k 1,2 ../commonSites/common4.txt > ../commonSites/sorted_common4.txt
time join -t# ../commonSites/sorted_common4.txt <(cat mcSites_HPNE_all.txt | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' | sed 's/ /\t/g') >| common4_mcSites_HPNE.txt
time join -t# ../commonSites/sorted_common4.txt <(cat mcSites_JD_all.txt | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' | sed 's/ /\t/g') >| common4_mcSites_JD.txt


#12-13-16
#now filter all the cpg sites
for chr in Y {22..1} X
do
    #mcSites are separated by spaces, not tabs
    echo $chr
    qsub ../filterCpg/filterCpg.sh mcSites_fixed/mcSites_JD_${chr}.txt mcSites_fixed/cpg/mcSites_JD_${chr}_cpg.txt
    qsub ../filterCpg/filterCpg.sh mcSites_fixed/mcSites_HPNE_${chr}.txt mcSites_fixed/cpg/mcSites_HPNE_${chr}_cpg.txt
done


#############
#scratch##########
join -t# ../commonSites/common4.txt <(cat mcSites_HPNE_all.txt | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/')


cat mcSites_HPNE_all.txt | sort -k 1,2 | sed 's/\(\w*\W*\w*\)\(.*\)/\1#\2/' | head -1000 > t1.txt
head t1.txt | sed 's/ /\t/g'


join -t# ../commonSites/sorted_common4.txt <(cat t1.txt | sed 's/ \t/g')

head -1000 ../commonSites/common4.txt
head -10 ../commonSites/common4.txt



