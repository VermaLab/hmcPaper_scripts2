#change this to ls *.1_val_1*
for fq1 in `ls ../trimmed/split/J432_AHHH5NBCXX_Lane1_ACAGTG_HPNE-1_BS.1_val_1*`;
do
    fq2=`echo $fq1 | sed 's/\(.*\)1_val_1\(.*\)/\12_val_2\2/'`
    pname=`echo $fq1 | sed 's/.*\/\(.*\).1_val_1.fq.\(p.*\)/\1_\2/'`
    #echo $fq1
    #echo $fq2
    #echo $pname
    echo "./makeBismarkBatch.sh $pname $fq1 $fq2 | qsub"
done
