library("VariantAnnotation")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene





#1.  Want to look at the methylation percentages of
#the jd and hpne samples at cpg sites.
#amit said the meth percentages looked a little high
#for the jd sample(over half had 100%)
#the mc histograms of the hmcSites should be high
#in the BS sample, and low in the OxBS sample


options(stringsAsFactors=F)


#load in the merged data
x = read.table("hmcSites_fixed2/cpg/hmcSites_BOTH_all_cpg_common4.txt")

colnames(x) = c("chr", "pos", 
	"jd.bs.meth", "jd.bs.tot", "jd.oxbs.meth", "jd.oxbs.tot", "jd.or", "jd.pval",
	"hpne.bs.meth", "hpne.bs.tot", "hpne.oxbs.meth", "hpne.oxbs.tot", "hpne.or", "hpne.pval"
)

head(x)
tail(x)
dim(x)
mean(x$jd.bs.tot)
mean(x$jd.oxbs.tot)
mean(x$hpne.bs.tot)
mean(x$hpne.oxbs.tot)

#for jd and hpne
#how many sites have at least 4 reads
#amongst those sites, what is avearge depth?





#try to get measures of overall methylation between samples
#what percentage of cpg sites have methylation less than 50%
ix.jd.lt50 = (x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"] <= .50)
ix.hpne.lt50 = (x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"] <= .50)
mean(x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"] <= .50)
mean(x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"] <= .50)


#what percentage of cpg sites have methylation greater than 50%
ix.jd.gt50 = (x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"] > .50)
ix.hpne.gt50 = (x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"] > .50)
mean(x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"] > .50)
mean(x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"] > .50)

#what percentage of cpg sites have total methylation
ix.jd.100 = (x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"] == 1)
ix.hpne.100 = (x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"] == 1)
mean(ix.jd.100)
mean(ix.hpne.100)

mean(x[,"jd.pval"] < 0.05)
mean(x[,"hpne.pval"] < 0.05)


#get indices of the signif hmc sites
jd.hmc.ix = which(x[,"jd.pval"] < 0.05)
hpne.hmc.ix = which(x[,"hpne.pval"] < 0.05)

length(jd.hmc.ix) / nrow(x)
length(hpne.hmc.ix) / nrow(x)
length(jd.hmc.ix)
length(hpne.hmc.ix)


#write the site lists
write.table(x[jd.hmc.ix,1:2], file="jd_hmc_p05.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(x[hpne.hmc.ix,1:2], file="hpne_hmc_p05.txt", col.names=F, row.names=F, quote=F, sep="\t")

write.table(x[ix.jd.gt50,1:2], file="jd_mc_gt50.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(x[ix.hpne.gt50,1:2], file="hpne_mc_gt50.txt", col.names=F, row.names=F, quote=F, sep="\t")

write.table(x[ix.jd.lt50,1:2], file="jd_mc_lt50.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(x[ix.hpne.lt50,1:2], file="hpne_mc_lt50.txt", col.names=F, row.names=F, quote=F, sep="\t")

write.table(x[ix.jd.100,1:2], file="jd_mc_100.txt", col.names=F, row.names=F, quote=F, sep="\t")
write.table(x[ix.hpne.100,1:2], file="hpne_mc_100.txt", col.names=F, row.names=F, quote=F, sep="\t")





#get percentage of sites that lie within exons introns promoters etc

getSiteCounts <- function(x){
	input = GRanges(seqnames=x[,1], ranges=IRanges(x[,2], x[,2]+1), strand="*")
	loc_hg19 <- locateVariants(input, txdb_hg19, AllVariants())
	loc_hg19 = loc_hg19[!duplicated(loc_hg19)]
	table(loc_hg19$LOCATION)
}

sCounts.jd.hmc = getSiteCounts(x[jd.hmc.ix,])
sCounts.hpne.hmc = getSiteCounts(x[hpne.hmc.ix,])

sCounts.jd.hmc
sCounts.hpne.hmc

sCounts.jd.hmc/sum(sCounts.jd.hmc)
sCounts.hpne.hmc/sum(sCounts.hpne.hmc)


#meth > 50% 
sCounts.jd.mc = getSiteCounts(x[ix.jd.gt50,])
sCounts.hpne.mc = getSiteCounts(x[ix.hpne.gt50,])
sCounts.jd.mc
sCounts.hpne.mc
sCounts.jd.mc/sum(sCounts.jd.mc)
sCounts.hpne.mc/sum(sCounts.hpne.mc)



#histograms of methylation across samples
png("methPerc_allCpg.png")
par(mfcol=c(2,2))
hist(x[,"jd.bs.meth"]/x[,"jd.bs.tot"], 100,         main="all cpg\njd bs meth %"    )
hist(x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"], 100,     main="all cpg\njd oxbs meth %"  )
hist(x[,"hpne.bs.meth"]/x[,"hpne.bs.tot"], 100,     main="all cpg\nhpne bs meth %"  )
hist(x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"], 100, main="all cpg\nhpne oxbs meth %")
dev.off()


png("methPerc_hmcCpg.png")
par(mfcol=c(2,2))
hist(x[jd.hmc.ix,"jd.bs.meth"]/x[jd.hmc.ix,"jd.bs.tot"], 100             ,  main="hmc cpg\njd bs meth %"    )
hist(x[jd.hmc.ix,"jd.oxbs.meth"]/x[jd.hmc.ix,"jd.oxbs.tot"], 100         ,  main="hmc cpg\njd oxbs meth %"  )
hist(x[hpne.hmc.ix,"hpne.bs.meth"]/x[hpne.hmc.ix,"hpne.bs.tot"], 100     ,  main="hmc cpg\nhpne bs meth %"  )
hist(x[hpne.hmc.ix,"hpne.oxbs.meth"]/x[hpne.hmc.ix,"hpne.oxbs.tot"],  100,  main="hmc cpg\nhpne oxbs meth %")
dev.off()


#get hydroxy meth percentages in and out of hmcsites
#hmc% is meth% bs - meth% oxbs

jd.hmcPerc = x[,"jd.bs.meth"]/x[,"jd.bs.tot"] - x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"]
jd.hmcPerc.hmcSites = x[jd.hmc.ix,"jd.bs.meth"]/x[jd.hmc.ix,"jd.bs.tot"] - x[jd.hmc.ix,"jd.oxbs.meth"]/x[jd.hmc.ix,"jd.oxbs.tot"]
hpne.hmcPerc = x[,"hpne.bs.meth"]/x[,"hpne.bs.tot"] - x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"]
hpne.hmcPerc.hmcSites = x[hpne.hmc.ix,"hpne.bs.meth"]/x[hpne.hmc.ix,"hpne.bs.tot"] - x[hpne.hmc.ix,"hpne.oxbs.meth"]/x[hpne.hmc.ix,"hpne.oxbs.tot"]

#how many sites with hmc > 50%?
mean(jd.hmcPerc > .50)
mean(hpne.hmcPerc > .50)
mean(x[,"jd.pval"] < 0.05)
mean(x[,"jd.pval"] < 0.01)

sum(jd.hmcPerc > 0 )
sum(jd.hmcPerc < 0 )
length(jd.hmc.ix)


jd.mcPerc.hmcSites

png("hmcPerc_cpg.png", 1000, 1000)
par(mfcol=c(2,2))
hist(jd.hmcPerc, 100)
hist(jd.hmcPerc.hmcSites, 100)
hist(hpne.hmcPerc, 100)
hist(hpne.hmcPerc.hmcSites, 100)
dev.off()


#overlay histograms of MC percentages
#with HMC percentages (of significant sites)

jd.mcPerc = x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"]
hpne.mcPerc = x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"]
jd.mcPerc.hmcSites = x[jd.hmc.ix,"jd.oxbs.meth"]/x[jd.hmc.ix,"jd.oxbs.tot"]
hpne.mcPerc.hmcSites = x[hpne.hmc.ix,"hpne.oxbs.meth"]/x[hpne.hmc.ix,"hpne.oxbs.tot"]


mean(jd.mcPerc)
sd(jd.mcPerc)
mean(hpne.mcPerc)
sd(hpne.mcPerc)

mean(jd.hmcPerc.hmcSites)
sd(jd.hmcPerc.hmcSites)
mean(hpne.hmcPerc.hmcSites)
sd(hpne.hmcPerc.hmcSites)


#examine the effect of smoothing on density plots
png("test1.png")
par(mfcol=c(3, 1))
plot(density(jd.mcPerc), col="red")
lines(density(jd.hmcPerc.hmcSites), col="blue")

plot(density(jd.mcPerc), col="red", ylim=c(0, 4))
lines(density(jd.hmcPerc.hmcSites), col="blue")

plot(density(jd.mcPerc, width=0.1),ylim=c(0, 4),  col="red")
lines(density(jd.hmcPerc.hmcSites, width=0.1), col="blue")
dev.off()



png("jd_histMCHmc_smooth.png")
plot(density(jd.mcPerc, width=0.1), ylim=c(0, 4), col="red", main="Methylation Percentage Densities:  JD", xlab="Methylation Percentage")
lines(density(jd.hmcPerc.hmcSites, width=0.1), col="blue")
legend("topleft", c("MC", "HMC"), pch=16, col=c("red", "blue"))
dev.off()

png("hpne_histMCHmc_smooth.png")
plot(density(hpne.mcPerc, width=0.1), ylim=c(0, 4), col="red", main="Methylation Percentage Densities:  HPNE", xlab="Methylation Percentage")
lines(density(hpne.hmcPerc.hmcSites, width=0.1), col="blue")
legend("topleft", c("MC", "HMC"), pch=16, col=c("red", "blue"))
dev.off()

png("hists_McPerc_smooth.png")
plot(density(jd.mcPerc, width=0.1), col="red", main="MC Percentage Densities", xlab="Methylation Percentage")
lines(density(hpne.mcPerc, width=.1), col="blue")
legend("topleft", c("JD", "HPNE"), pch=16, col=c("red", "blue"))
dev.off()

png("hists_HmcPerc_smooth.png")
plot(density(jd.hmcPerc.hmcSites, width=0.10), col="red", ylim=c(0, 4), main="HMC Percentage Densities", xlab="Hydroxy-Methylation Percentage")
lines(density(hpne.hmcPerc.hmcSites, width=.10), col="blue")
legend("topleft", c("JD", "HPNE"), pch=16, col=c("red", "blue"))
dev.off()


#JD
png("jd_histMcHmc.png")
#methylation histogram over all cpg sites
hist.mc = hist(x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"], plot=F , 100 )
#hmc histogram over all significant cpg sites
hist.hmc.sig = hist(jd.hmcPerc.hmcSites, plot=F, 100)
xlim <- range(hist.mc$breaks,hist.hmc.sig$breaks)
ylim <- range(0,hist.mc$density, hist.hmc.sig$density)
## plot the first graph
plot(hist.mc,xlim = xlim, ylim = ylim,
     col = rgb(1,0,0,0.4),xlab = 'Lengths',
     freq = FALSE, ## relative, not absolute frequency
     main = 'JD:  Methylation')
## plot the second graph on top of this
opar <- par(new = FALSE)
plot(hist.hmc.sig,xlim = xlim, ylim = ylim,
     xaxt = 'n', yaxt = 'n', ## don't add axes
     col = rgb(0,0,1,0.4), add = TRUE,
     freq = FALSE) ## relative, not absolute frequency
## add a legend in the corner
legend('topleft',c('MC','HMC'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
par(opar)
dev.off()



#HPNE
png("hpne_histMcHmc.png")
#methylation histogram over all cpg sites
hist.mc = hist(x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"], plot=F , 100 )
#hmc histogram over all significant cpg sites
hist.hmc.sig = hist(hpne.hmcPerc.hmcSites, plot=F, 100)
xlim <- range(hist.mc$breaks,hist.hmc.sig$breaks)
ylim <- range(0,hist.mc$density, hist.hmc.sig$density)
## plot the first graph
plot(hist.mc,xlim = xlim, ylim = ylim,
     col = rgb(1,0,0,0.4),xlab = 'Lengths',
     freq = FALSE, ## relative, not absolute frequency
     main = 'HPNE:  Methylation')
## plot the second graph on top of this
opar <- par(new = FALSE)
plot(hist.hmc.sig,xlim = xlim, ylim = ylim,
     xaxt = 'n', yaxt = 'n', ## don't add axes
     col = rgb(0,0,1,0.4), add = TRUE,
     freq = FALSE) ## relative, not absolute frequency
## add a legend in the corner
legend('topleft',c('MC','HMC'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
par(opar)
dev.off()



#overlay hmc jd and hpne
png("hists_HmcPerc.png")
hist.hmc.jd= hist(jd.hmcPerc.hmcSites, plot=F, 100)
hist.hmc.hpne= hist(hpne.hmcPerc.hmcSites, plot=F, 100)
xlim <- range(hist.hmc.jd$breaks,hist.hmc.hpne$breaks)
ylim <- range(0,hist.hmc.jd$density, hist.hmc.hpne$density)
## plot the first graph
plot(hist.hmc.jd,xlim = xlim, ylim = ylim,
     col = rgb(1,0,0,0.4),xlab = 'Lengths',
     freq = FALSE, ## relative, not absolute frequency
     main = 'HMC: JD, HPNE')
## plot the second graph on top of this
opar <- par(new = FALSE)
plot(hist.hmc.hpne,xlim = xlim, ylim = ylim,
     xaxt = 'n', yaxt = 'n', ## don't add axes
     col = rgb(0,0,1,0.4), add = TRUE,
     freq = FALSE) ## relative, not absolute frequency
## add a legend in the corner
legend('topleft',c('JD','HPNE'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
par(opar)
dev.off()





#overlay mc: jd, hpne
png("hists_McPerc.png")
hist.mc.jd = hist(x[,"jd.oxbs.meth"]/x[,"jd.oxbs.tot"], plot=F, 100  )
hist.mc.hpne = hist(x[,"hpne.oxbs.meth"]/x[,"hpne.oxbs.tot"], plot=F, 100)
xlim <- range(hist.mc.jd$breaks,hist.mc.hpne$breaks)
ylim <- range(0,hist.mc.jd$density, hist.mc.hpne$density)
## plot the first graph
plot(hist.mc.jd,xlim = xlim, ylim = ylim,
     col = rgb(1,0,0,0.4),xlab = 'Lengths',
     freq = FALSE, ## relative, not absolute frequency
     main = 'MC: JD, HPNE')
## plot the second graph on top of this
opar <- par(new = FALSE)
plot(hist.mc.hpne,xlim = xlim, ylim = ylim,
     xaxt = 'n', yaxt = 'n', ## don't add axes
     col = rgb(0,0,1,0.4), add = TRUE,
     freq = FALSE) ## relative, not absolute frequency
## add a legend in the corner
legend('topleft',c('JD','HPNE'),
       fill = rgb(1:0,0,0:1,0.4), bty = 'n',
       border = NA)
par(opar)
dev.off()






# there are a total of 27,999,538 cpg sites in genome
# we have data for 11,296,328 of them



#brd promoter region
#look at hmc on chr19 between these points
p1 = 15391262 - 1500
p2 = 15443342 + 1500
head(x)
x.19 = x[x[,1] == "chr19",]
ix.brd4 = which((x.19$pos >= p1) & (x.19$pos < p2))
length(ix.brd4)
#plot out mc and hmc over this region
x.brd4 = x.19[ix.brd4,]
head(x.brd4)
jd.mc = x.brd4$jd.oxbs.meth/x.brd4$jd.oxbs.tot
jd.hmc = x.brd4[,"jd.bs.meth"]/x.brd4[,"jd.bs.tot"] - x.brd4[,"jd.oxbs.meth"]/x.brd4[,"jd.oxbs.tot"]
hpne.mc = x.brd4$hpne.oxbs.meth/x.brd4$hpne.oxbs.tot
hpne.hmc = x.brd4[,"hpne.bs.meth"]/x.brd4[,"hpne.bs.tot"] - x.brd4[,"hpne.oxbs.meth"]/x.brd4[,"hpne.oxbs.tot"]
hmc
stem(hmc)

#draw the less confidence dots with transparency
blues = sapply(1:100, function(i){
	rgb(0,0,1,i/100.0)
})
jd.hmc.cols =blues[101 - round(x.brd4$jd.pval, 2)*100]
hpne.hmc.cols =blues[101 - round(x.brd4$hpne.pval, 2)*100]


reds = sapply(1:100, function(i){
	rgb(1,0,0,i/100.0)
})
#more reads -> more confidences
n = max(x.brd4$jd.oxbs.tot)
jd.mc.cols = reds[round(x.brd4$jd.oxbs.tot / n * 100)]
hpne.mc.cols = reds[round(x.brd4$hpne.oxbs.tot / n * 100)]

png("jd_brd4_mc_hmc_conf.png")
plot(x.brd4$pos, jd.hmc, col=jd.hmc.cols, ylim=c(-1, 1), ylab="methylation")
points(x.brd4$pos, jd.mc, col=jd.mc.cols)
dev.off()

png("hpne_brd4_mc_hmc_conf.png")
plot(x.brd4$pos, hpne.hmc, col=hpne.hmc.cols, ylim=c(-1, 1), ylab="methylation")
points(x.brd4$pos, hpne.mc, col=hpne.mc.cols)
dev.off()

png("jd_brd4_mc_hmc.png")
plot(x.brd4$pos, hmc, col="blue", ylim=c(-1, 1), ylab="methylation")
points(x.brd4$pos, mc, col="red")
dev.off()

x11()
plot(x.brd4$pos, hpne.hmc, col="blue", ylim=c(-1, 1), ylab="methylation")
points(x.brd4$pos, hpne.mc, col="red")



plot(x.brd4$pos, ylim=c(-1, 1), xlim=c(min(x.brd4$pos), max(x.brd4$pos)))
segments(x0=x.brd4$pos, y0=x.brd4$pos*0, y1=mc, col="red")
segments(x0=x.brd4$pos, y0=x.brd4$pos*0, y1=hmc, col="blue")

#find places where methylation is significantly different between
#hpne and jd
diffmeth.pvals = NA*numeric(nrow(x))
diffmeth.or = NA*numeric(nrow(x))
for (i in 1:nrow(x)){
	if (i %% 10000 == 0){
		print(i)
	}

	#a = 10
	#b = 20
	#c = 30
	#d = 40
	a = x[i,"jd.bs.meth"] #num converted reads
	b = x[i,"jd.bs.tot"] - a #num of converted reads
	c = x[i,"hpne.bs.meth"]
	d = x[i,"hpne.bs.tot"] - c

	#see if there's a significant different in the ratios
	#of methylated reads between jd and hpne
	res = fisher.test(matrix(c(a, b ,c, d), ncol=2))
	p = res$p.value
	est = res$estimate
	meth.pvals[i] = p
	meth.or[i] =est 
}

save.image("work_diffMeth.rData")
qvals = p.adjust(meth.pvals)

ix.sig =  qvals < 0.05

#is this hmc loss and gain?
#more meth in JD, sig hmc in HPNE
sum(ix.sig & meth.or > 1 & x[,"hpne.pval"] < 0.05)
head(x[ix.sig & meth.or > 1 & x[,"hpne.pval"] < 0.05,])

#more meth in hpne, sig hmc in JD
sum(ix.sig & meth.or < 1 & x[,"jd.pval"] < 0.05)
x[ix.sig & meth.or < 1 & x[,"jd.pval"] < 0.05,]








#collect the four sets
#sites with significant HMC in JD
#sites with significant HMC in HPNE 
#sites with significantly more MC in JD than HPNE
#sites with significantly more MC in HPNE than JD



#get site distributions for these four sets
sum(x[,"jd.pval"] < 0.05)
sum(x[,"jd.pval"] < 0.01)
sum(x[,"hpne.pval"] < 0.01)
sum(x[,"hpne.pval"] < 0.05)

png("methSigDiff_hist.png")
par(mfcol=c(2,1))
hist(log(meth.or), 1000)
hist(log(meth.or[qvals < 0.2]), 1000)
dev.off()

qvals = p.adjust(meth.pvals)
sum(meth.or < 1 & qvals < 0.2)
sum(meth.or > 1 & qvals < 0.2)
sum(qvals < 0.2)
ix = meth.or > 1 & qvals < 0.2
head(x[ix,])


jd.mc.sites = x[meth.or > 1 & qvals < 0.2,1:2]
hpne.mc.sites = x[meth.or < 1 & qvals < 0.2,1:2]
jd.hmc.sites = x[(x[,"jd.or"] > 1) & (x[,"jd.pval"] < 0.05),1:2]
hpne.hmc.sites = x[(x[,"hpne.or"] > 1) & (x[,"hpne.pval"] < 0.05),1:2]
dim(jd.mc.sites)
dim(hpne.mc.sites)
dim(jd.hmc.sites)
dim(hpne.hmc.sites)


write.table(file="sites_fixed2_jdMc.txt", jd.mc.sites, col.names=F, row.names=F, quote=F, sep="\t")
write.table(file="sites_fixed2_hpneMc.txt", hpne.mc.sites, col.names=F, row.names=F, quote=F, sep="\t")
write.table(file="sites_fixed2_jdHmc_p05.txt", jd.hmc.sites, col.names=F, row.names=F, quote=F, sep="\t")
write.table(file="sites_fixed2_hpneHmc_p05.txt", hpne.hmc.sites, col.names=F, row.names=F, quote=F, sep="\t")


#what if we only look at the sites with significant meth
x[,"jd.meth"]

i = 1

a = x[i,"jd.bs.meth"] 
b = x[i,"jd.bs.tot"] - a
c = x[i,"hpne.bs.meth"] 
d = x[i,"hpne.bs.tot"] - c

matrix(c(a,b,c,d), 2)

ix = which((qvals < 0.05) & (x[,"jd.pval"] < 0.05))
a = x[ix, "jd.pval"]
b = p.adjust(a)
c = which(b < 0.05)
x[ix[c],]
which(p.adjust(x[(meth.pvals < 0.05) & (x[,"jd.pval"] < 0.05), "jd.pval"]) < 0.05)
x[c(19880, 24423),]
sum((meth.pvals > 0.05) & (x[,"jd.pval"] < 0.05))
sum((qvals >= 0.05) & (x[,"jd.pval"] < 0.05))

sum((qvals < 0.05) & (x[qvals < 0.05,"hpne.pval"] < 0.05))
q1 = p.adjust(x[qvals < 0.05,"hpne.pval"])
sum(q1 < 0.15)

sum(meth.pvals < 0.05)
fisher.test(








######################################
##scratch
#view the histograms of odds ratios.

#write out the places with significant difference between jd and hpne

n = 320000
hist(meth.pvals, 100)
hist(meth.pvals, 100)

hist(log(meth.or[meth.or != 1]), 1000)
hist(log(meth.or[meth.or != 1]), 100)

hist(log(meth.or[1:n][meth.or[1:n] != 1]), 1000)
hist(log(meth.or[1:n][meth.or[1:n] != 1]), 100)
sum(meth.pvals[1:n] < 0.05, na.rm=T)
qvals = p.adjust(meth.pvals)
mean(qvals < 0.05)
sum(qvals < 0.05)
hist(log(meth.or[qvals < 0.05]), 100)
ix = which(qvals < 0.05)
length(ix)
x[ix,]
