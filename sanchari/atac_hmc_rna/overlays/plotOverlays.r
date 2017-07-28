#1.  jd atac hmc profile
#       overlayed with
#     jd atac no hmc profile
#
#
#2.  hpne atac hmc profile
#       overlayed with
#     hpne atac no hmc profile
#



#plot(regcovMat, type='l', xaxt="n")
#axis(side=1, at=xticks$pos, labels = xticks$lab)


#load in the jd peak sets

load("jd_strandedAtacPeaks_hmc/avgprof.RData")
jd.hmc.regcovMat = regcovMat
jd.hmc.xticks = xticks

load("jd_strandedAtacPeaks_nohmc/avgprof.RData")
jd.nohmc.regcovMat = regcovMat
jd.nohmc.xticks = xticks


png("jd_atac_rna_hmc.png")
plot(jd.hmc.regcovMat, type='l', xaxt="n", xlab="Genomic Region (5' −> 3')", ylab="Read count Per Million mapped reads", col="blue", lwd=3)
lines(jd.nohmc.regcovMat, type='l', xaxt="n", col="red", lw=3)
axis(side=1, at=jd.hmc.xticks$pos+1, labels = jd.hmc.xticks$lab)
abline(v=floor(pts/2)+1, col=rgb(.5,.5,.5))
legend("topleft", c("jd.hmc.atac.rna", "jd.nohmc.atac.rna"), col=c("blue", "red"), pch=16)
dev.off()


load("hpne_strandedAtacPeaks_hmc/avgprof.RData")
hpne.hmc.regcovMat = regcovMat
hpne.hmc.xticks = xticks

load("hpne_strandedAtacPeaks_nohmc/avgprof.RData")
hpne.nohmc.regcovMat = regcovMat
hpne.nohmc.xticks = xticks


png("hpne_atac_rna_hmc.png")
plot(hpne.hmc.regcovMat, type='l', xaxt="n", xlab="Genomic Region (5' −> 3')", ylab="Read count Per Million mapped reads", col="blue", lwd=3)
lines(hpne.nohmc.regcovMat, type='l', xaxt="n", col="red", lw=3)
axis(side=1, at=hpne.hmc.xticks$pos+1, labels = hpne.hmc.xticks$lab)
abline(v=floor(pts/2)+1, col=rgb(.5,.5,.5))
legend("topleft", c("hpne.hmc.atac.rna", "hpne.nohmc.atac.rna"), col=c("blue", "red"), pch=16)
dev.off()
