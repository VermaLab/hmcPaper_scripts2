install.packages("Gviz")
library("Gviz")
install.packages("playwith")
library("playwith")
options(stringsAsFactors=F)


region=data.frame(name = "brd4", chr="19", start=15346301, end=15393262)
region$start = region$start - 2000
region$end  = region$end + 2000


writeGeneInfo <- function(name, chr, start, end, W = 2000){
	region=data.frame(name = name, chr=chr, start=start, end=end)
	region$start = region$start - W
	region$end  = region$end + W

	#where the conversion counts are stored
	ccfol = "/home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/conversionCounts/"

	#load in the jd and hpne info for the region
	#takes about 3 minutes to load in chr19
	t1 = proc.time()
	jd.oxbs = read.table(paste0(ccfol, "JD_OxBs_", region$chr, ".txt.gz"))
	jd.bs = read.table(paste0(ccfol, "JD_BS_", region$chr, ".txt.gz"))
	hpne.oxbs = read.table(paste0(ccfol, "HPNE_OxBS_", region$chr, ".txt.gz"))
	hpne.bs = read.table(paste0(ccfol, "HPNE_BS_", region$chr, ".txt.gz"))
	t2 = proc.time()
	print(t2 - t1)

	#get the rows that lie within the genes region
	jd.oxbs.reg = jd.oxbs[head(which(jd.oxbs[,2] > region$start),1):head(which(jd.oxbs[,2] > region$end),1),]
	jd.bs.reg = jd.bs[head(which(jd.bs[,2] > region$start),1):head(which(jd.bs[,2] > region$end),1),]
	hpne.oxbs.reg = hpne.oxbs[head(which(hpne.oxbs[,2] > region$start),1):head(which(hpne.oxbs[,2] > region$end),1),]
	hpne.bs.reg = hpne.bs[head(which(hpne.bs[,2] > region$start),1):head(which(hpne.bs[,2] > region$end),1),]


	#
	#merge together bs and oxbs
	jd.reg = merge(jd.bs.reg, jd.oxbs.reg, by=2)
	hpne.reg = merge(hpne.bs.reg, hpne.oxbs.reg, by=2)

	#playwith({
	#	plot(jd.reg[,1], jd.reg[,3]/jd.reg[,4], type="l")
	#})

	write.table(paste0(region$name, "_jd_bs.bdg"), x=(cbind(jd.reg[,2], jd.reg[, 1], jd.reg[,1]+1, jd.reg[,3]/jd.reg[,4])), sep="\t", col.names=F, row.names=F, quote=F)
	write.table(paste0(region$name, "_jd_oxbs.bdg"), x=(cbind(jd.reg[,2], jd.reg[, 1], jd.reg[,1]+1, jd.reg[,6]/jd.reg[,7])), sep="\t", col.names=F, row.names=F, quote=F)

	write.table(paste0(region$name, "_hpne_bs.bdg"), x=(cbind(hpne.reg[,2], hpne.reg[, 1], hpne.reg[,1]+1, hpne.reg[,3]/hpne.reg[,4])), sep="\t", col.names=F, row.names=F, quote=F)
	write.table(paste0(region$name, "_hpne_oxbs.bdg"), x=(cbind(hpne.reg[,2], hpne.reg[, 1], hpne.reg[,1]+1, hpne.reg[,6]/hpne.reg[,7])), sep="\t", col.names=F, row.names=F, quote=F)

}




#brd4: chr19:15,346,301-15,393,262
#myc:  chr8:128,746,315-128,755,680
#vegfa: chr6:43,735,946-43,756,223
#vegfb: chr11:64,000,056-64,008,736
#vegfc: chr4:177,602,689-177,715,899

writeGeneInfo("brd4", "19", 15346301, 15393262)
writeGeneInfo("myc", "8", 128746315, 128755680)
writeGeneInfo("vegfa", "6", 43735946, 43756223)
writeGeneInfo("vegfb", "11", 64000056, 64008736)
writeGeneInfo("vegfc", "4", 177602689, 177715899)



