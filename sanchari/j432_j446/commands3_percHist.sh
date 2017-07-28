#find the percetages of hmc at each selected site

cat hmcSites/hmcSites_JD_all.txt | awk '{print $5/$6 - $3/$4}' > hmcPerc_JD.txt
cat hmcSites/hmcSites_HPNE_all.txt | awk '{print $5/$6 - $3/$4}' > hmcPerc_HPNE.txt


#R code
x = read.table("hmcPerc_JD.txt")
x = x[,1]
png("hmcPerc_JD.png")
hist(x, 100, main="hmC Percentage of Selected Sites: JD")
dev.off()

x = read.table("hmcPerc_HPNE.txt")
x = x[,1]
png("hmcPerc_HPNE.png")
hist(x, 100, main="hmC Percentage of Selected Sites: HPNE")
dev.off()



join <(zcat conversionCounts/HPNE_BS_22.txt.gz | cut -f2-4) <(zcat conversionCounts/HPNE_OxBS_22.txt.gz | cut -f2-4) \
 | awk '{if ($3 >= 4 && $5 >= 4) print $4/$5 - $2/$3}' > commonPerc_HPNE22.txt


#r code
x = read.table("commonPerc_HPNE22.txt")
x = x[,1]
png("commonPerc_HPNE22.png")
hist(x, 100, main="hmC Percentage of Common Sites: HPNE chr22")
dev.off()





