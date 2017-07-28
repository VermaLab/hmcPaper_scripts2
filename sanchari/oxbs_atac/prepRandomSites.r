#read in a bed file
options(stringsAsFactors=F)

#sites = read.table("sites_hpne_all.bed", header=F)
sites = read.table("commonSites1_hpne.bed", header=F)

dim(sites)
head(sites)

#for each chrom
sites.all = lapply(unique(sites[,1]), function(chrom){
    ix = sites[,1] == chrom
    sites[ix,]
})
names(sites.all) = unique(sites[,1])
length(sites.all)

#for each chromosome write random sites between min and max
for (chrom in unique(sites[,1])){
    print(chrom)
    my.sites = sites.all[[chrom]]

    r1 = min(my.sites[,2])
    r2 = max(my.sites[,2])

    #range(min/max)
    rmy.sites = my.sites
    rmy.sites[,2] = sample(r1:r2, nrow(rmy.sites))
    rmy.sites[,3] = rmy.sites[,2]+1

    
    write.table(file="sites_rand_all.bed", x=rmy.sites, row.names=F, col.names=F, quote=F, sep="\t", append=T)
}




#replace coordinates with random numbers in the same
r1 = min(sites.19[,2])
r2 = max(sites.19[,2])

#range(min/max)
rsites.19 = sites.19
rsites.19[,2] = sample(r1:r2, nrow(rsites.19))
rsites.19[,3] = rsites.19[,2]+1

#write new bed file
write.table(file="sites_rand_chr19.bed", x=rsites.19, row.names=F, col.names=F, quote=F, sep="\t")



#prep sites that appear in one sample but not another.

sites.hpne = read.table("sites_hpne_all.bed", header=F)
sites.jd = read.table("sites_jd_all.bed", header=F)

sapply(paste0("chr",c(1:22,"X","Y")), function(chrom){
	in1 = (sites.hpne[sites.hpne[,1] == chrom,2] %in% sites.jd[sites.jd[,1] == chrom,2])
})



