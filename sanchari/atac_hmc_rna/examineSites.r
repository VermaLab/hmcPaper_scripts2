options(stringsAsFactors=F)

#look at the atac-seq sites
#	separate into 2 sets, those that have hmc sites within
#	those that don't
#	maybe break into sets based on the percentage of hmc sites

#compare these two site sites to rna expression
#ngs.plot


#will need to add strand info



#load in the atac seq peaks
atac.jd = read.table("jd_strandedAtacPeakSummits.txt")
atac.hpne = read.table("hpne_strandedAtacPeakSummits.txt")
head(atac.jd)
head(atac.hpne)


#load in the hmc sites
hmc.jd = read.table("sites_fixed2_jdHmc_p05.txt")
hmc.hpne = read.table("sites_fixed2_hpneHmc_p05.txt")
head(hmc.jd)
head(hmc.hpne)


#how many sites have hmc sites within them?
x = atac.jd
y = hmc.jd

#work chrom by chrom
w = 3000
unique(x[,1])

#x are the atac peak summits
#y are the hmc sites
#w is the width
#for each atac peak, return the number of hmc sites within w bases 
getOverlapPeaks <- function(x, y, w =3000){
	peaks = sapply(1:nrow(x), function(i){
		if (i %% 1000 == 0){
			print (i)
		}
		ix = y[,1] == x[i,1]
		sum((y[ix,2] > x[i,2]-w) & (y[ix,2] < x[i,3]+w))
	})
	peaks
}

hmcPeaks.jd = getOverlapPeaks(atac.jd, hmc.jd, 3000)
hmcPeaks.hpne = getOverlapPeaks(atac.hpne, hmc.hpne, 3000)

write.table(file="jd_strandedAtacPeaks_hmc.txt", x=atac.jd[hmcPeaks.jd > 0,], row.names=F, col.names=F, quote=F, sep="\t")
write.table(file="jd_strandedAtacPeaks_nohmc.txt", x=atac.jd[hmcPeaks.jd == 0,], row.names=F, col.names=F, quote=F, sep="\t")
write.table(file="hpne_strandedAtacPeaks_hmc.txt", x=atac.hpne[hmcPeaks.hpne > 0,], row.names=F, col.names=F, quote=F, sep="\t")
write.table(file="hpne_strandedAtacPeaks_nohmc.txt", x=atac.hpne[hmcPeaks.hpne == 0,], row.names=F, col.names=F, quote=F, sep="\t")

dim(x)
dim(y)


#separate the atac sites into those with/without hmc sites

