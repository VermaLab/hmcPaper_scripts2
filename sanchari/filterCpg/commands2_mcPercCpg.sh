#get the average methylation percentage
#regular C(not CpG) sites
#CpG sites



#run the bedgraph output through cpg filter
fol=~/mnt/hpc_home/projects/sanchari/j432_j446/conversionCounts/
ls $fol


file=${fol}JD_BS_22.txt.gz
echo $file


#separate the C's from CpG's
for file in `ls ${fol}*JD_BS*.gz`
do
    echo $file
    #The CpG sites
    echo "cpg filtering..."
    zcat $file | python ../filterCpg.py | awk '{printf "%s\t%d#%d\t%d\n", $1,$2,$3,$4}' >| temp_cpg.txt
    #all sites
    zcat $file | awk '{printf "%s\t%d#%d\t%d\n", $1,$2,$3,$4}' >| temp_allC.txt
    #append just the C sites
    echo "writing C's..."
    join -v1 -t"#" <(sort -k1,1 -t"#" temp_allC.txt) <(sort -k1,1 -t"#" temp_cpg.txt) >> JD_BS_C.txt
    #append the CpG sites
    echo "writing CpG's..."
    cat temp_cpg.txt >> JD_BS_CpG.txt
done

#now compute average methylation between C's and CpG's


cat Cpgfiles/JD_BS_CpG.txt | ./calcRunningAve.pl >| meanCpg.txt
cat Cpgfiles/JD_BS_C.txt | ./calcRunningAve.pl >| meanC.txt 

