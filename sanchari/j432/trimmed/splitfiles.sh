#split each file into 40M line chunks


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



