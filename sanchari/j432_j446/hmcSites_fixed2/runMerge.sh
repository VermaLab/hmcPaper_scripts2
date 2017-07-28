
for chr in Y {22..1} X
do
    echo $chr
    Rscript mergeSites.r $chr
done
