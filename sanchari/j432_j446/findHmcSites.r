#!/usr/bin/env Rscript

options(stringsAsFactors=F)

#command line arguments
args = commandArgs(trailingOnly=TRUE)
f.bs = args[1]
f.oxbs = args[2]
f.hmc = args[3]
#f.bs = "conversionCounts/HPNE_BS_22.txt.gz" 
#f.oxbs = "conversionCounts/HPNE_OxBS_22.txt.gz"
#f.hmc = "hmcSites_HPNE_22.txt"

print(f.bs)
print(f.oxbs)
print(f.hmc)

#BS
#c -> t
#5mc -> c
#5hmc -> c
#
#oxBS
#c -> t
#5mc -> c
#5hmc -> t

#we want to find the sites with hmc
#low conversion in BS
#high conversion in oxBS


#candidate sites should have 
#conversion rate in Bs  below 50%
#conversion rate in oxBs above 50%



#read in the count info from the gz files
print(paste0("loading bs:  ", f.bs))
x1 = read.table(pipe(paste0("zcat ", f.bs)))
print(paste0("loading oxbs:  ", f.oxbs))
x2 = read.table(pipe(paste0("zcat ", f.oxbs)))

#head(x1)
#merge the bs and oxbs together
print("merging bs/oxbs")
x = merge(x1, x2, by.x=2, by.y=2)


#for each position calculate probability that
#oxbs conversion is more than bs conversion
#these are the hmc sites

#store data like this will make access faster
vals = t(data.matrix(x[,c("V3.x", "V4.x", "V3.y", "V4.y")]))

print ("testing...")
t1 = Sys.time()
#res.sites = sapply(1:20000, function(i){
res.sites = sapply(1:nrow(x), function(i){
	if (i %% 1000 == 0){
		print(i)
	}
	v = vals[,i]
	#v = unlist(x[i,c("V3.x", "V4.x", "V3.y", "V4.y")])
	m = matrix((c(v[1], v[3], v[2]-v[1], v[4]-v[3])), ncol=2)

	res = fisher.test(m)
	c(odds=res$estimate,pval=res$p.val)
})
t2 = Sys.time()
t2 - t1
dim(res.sites)

#lets just look at the sites with pval < 0.05
#and where oxbs has more conversion than bs
ix = which((res.sites[1,] < 1) & (res.sites[2,] < 0.05))
res = cbind(x[ix, 2:1], x[ix,c("V3.x", "V4.x", "V3.y", "V4.y")], t(res.sites[,ix]))
print (paste0("writing:  ", f.hmc))
write.table(x=res, file=f.hmc, row.names=F, sep="\t")



