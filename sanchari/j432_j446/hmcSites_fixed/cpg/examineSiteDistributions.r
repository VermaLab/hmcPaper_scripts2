library("VariantAnnotation")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

options(stringsAsFactors=F)



x = read.table("hmcSites_BOTH_all_cpg_common4.txt")


#look at the hmcsites with pval < 0.05
head(x)
png("jd_methylation_hist.png")
hist(x[,3]/x[,4], 100, xlab="methylation percentage")
dev.off()
mean(x[,3]/x[,4] > .50)
mean(x[,3]/x[,4] == 1)
dim(x)
#take the top 10000 with the most methylatiojh

ix = order(x[,3]/x[,4], x[,4])
#ix.mc2 = tail(ix, 10000)
ix.mc = tail(ix, 100000)
head(x[ix,])
ix.mc = sort(sample(which(x[,3]/x[,4] > .50), 10000 ))

head(x[ix.mc,])
tail(x[ix.mc,])
y = x[ix.mc,1:2]
y[,3] = "+"
write.table(y, file="test1_jdmc_top100k.txt", col.names=F, row.names=F, sep="\t", quote=F)

jd.hmc.ix = which(x[,8] < 0.01)
hpne.hmc.ix = which(x[,14] < 0.01)
length(jd.hmc.ix)
length(hpne.hmc.ix)

ov.hmc.ix = jd.hmc.ix %in% hpne.hmc.ix

sum(ov.hmc.ix)
head(x[jd.hmc.ix,])

y = x[jd.hmc.ix,1:2]
y[,3] = "+"
write.table(y, file="test1_jdhmc.txt", col.names=F, row.names=F, sep="\t", quote=F)

#write a set of random sites

jd.bad.ix = which(x[,8] > 0.91)
length(jd.bad.ix)
jd.bad.ix = sample(jd.bad.ix, length(jd.hmc.ix))
length(jd.bad.ix)
y = x[jd.bad.ix,1:2]
y[,3] = "+"
write.table(y, file="test1_jdbad.txt", col.names=F, row.names=F, sep="\t", quote=F)




#look at sites with high methylation





getSiteCounts <- function(x){
    input = GRanges(seqnames=x[,1], ranges=IRanges(x[,2], x[,2]+1), strand="*")
    loc_hg19 <- locateVariants(input, txdb_hg19, AllVariants())
    loc_hg19 = loc_hg19[!duplicated(loc_hg19)]
    table(loc_hg19$LOCATION)
}


sCounts.1 = getSiteCounts(x[jd.hmc.ix,1:2])
sCounts.1a = getSiteCounts(x[jd.bad.ix,1:2])
sCounts.2 = getSiteCounts(x[hpne.hmc.ix,1:2])

barplot(sCounts.1)
x11(); barplot(sCounts.1a)
x11(); barplot(sCounts.2)


