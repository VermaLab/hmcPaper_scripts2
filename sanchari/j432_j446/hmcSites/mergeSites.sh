#add the header to the merged file
head -1 hmcSites_HPNE_22.txt >| hmcSites_HPNE_all.txt
head -1 hmcSites_JD_22.txt >| hmcSites_JD_all.txt
for chr in {1..22} X Y
do
    echo $chr
    #process the HPNE files
    echo `ls *hmcSites_HPNE_${chr}.txt`
    file=`ls *hmcSites_HPNE_${chr}.txt`
    if echo $file | grep -e "stream"; then
        cat $file >> hmcSites_HPNE_all.txt
    else
        echo "nostream"
        #skip the first line and take out quotations
        tail -n +2 $file  | sed 's/"//g' >> hmcSites_HPNE_all.txt
    fi
    #process the JD files
    echo `ls *hmcSites_JD_${chr}.txt`
    file=`ls *hmcSites_JD_${chr}.txt`
    if echo $file | grep -e "stream"; then
        cat $file >> hmcSites_JD_all.txt
    else
        echo "nostream"
        #skip the first line and take out quotations
        tail -n +2 $file  | sed 's/"//g' >> hmcSites_JD_all.txt
    fi
done

