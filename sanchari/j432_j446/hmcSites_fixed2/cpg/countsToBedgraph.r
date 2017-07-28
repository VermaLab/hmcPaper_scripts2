options(stringsAsFactors=F)

x = read.table("hmcSites_BOTH_all_cpg_common4.txt")
head(x)


colnames(x) = c("chr", "pos", 
	"jd.bs.meth", "jd.bs.tot", "jd.oxbs.meth", "jd.oxbs.tot", "jd.or", "jd.pval",
	"hpne.bs.meth", "hpne.bs.tot", "hpne.oxbs.meth", "hpne.oxbs.tot", "hpne.or", "hpne.pval"
)

head(x)


#map the methylated reads to positive(blue)
#unmethylation reads to negative(red)

#do this for BS and OxBS



#positive methylated blue

descr="jd.bs"
col1="jd.bs.meth"
col2="jd.bs.tot"
#show raw counts of converted and non-converted reads
writeBedgraph_v1 <- function(descr, x, col1, col2){
	#turn off scientific notation
	options(scipen=999)

	ofile = paste0(descr, "_meth.bedGraph")
	trackinfo = paste0("track type=bedGraph name=\"",ofile,"\" visibility=full color=0,0,255 altColor=255,255,255")
	x.meth = cbind(x[,"chr"], x[,"pos"], x[,"pos"]+1, x[,col1])
	write(trackinfo, file = ofile)
	write.table(file= ofile, x.meth, col.names=F, row.names=F, quote=F, sep="\t", append=T)

	#negative non-methylation red
	ofile = paste0(descr, "_nometh.bedGraph")
	trackinfo = paste0("track type=bedGraph name=\"",ofile,"\" visibility=full altColor=255,0,0 color=255,255,255")
	x.nometh = cbind(x[,"chr"], x[,"pos"], x[,"pos"]+1, -(x[,col2] - x[,col1]))
	write(trackinfo, file = ofile)
	write.table(file= ofile, x.nometh, col.names=F, row.names=F, quote=F, sep="\t", append=T)

}


writeBedgraph_v1("jd.bs", x, "jd.bs.meth", "jd.bs.tot")
writeBedgraph_v1("jd.oxbs", x, "jd.oxbs.meth", "jd.oxbs.tot")
writeBedgraph_v1("hpne.bs", x, "hpne.bs.meth", "hpne.bs.tot")
writeBedgraph_v1("hpne.oxbs", x, "hpne.oxbs.meth", "hpne.oxbs.tot")





#now show the hmc and mc estimates as a bedgraph
samp = "jd"
writeBedgraph_v2 <- function(samp, x){
	#separate mc/hmc into two graphs

	#ratio of meth reads to total reads is estiate of methylation
	r.bs = x[,paste0(samp, ".bs.meth")]/x[,paste0(samp, ".bs.tot")]
	r.oxbs = x[,paste0(samp, ".oxbs.meth")]/x[,paste0(samp, ".oxbs.tot")]
	
	#hist(r.bs, 100)
	#hist(r.oxbs - r.bs, 100)

	#difference of BS ratio and OxBS ratio is estimate for hmc%
	#shown in blue
	ofile = paste0(samp, "_hmc.bedGraph")
	trackinfo = paste0("track type=bedGraph name=\"",ofile,"\" visibility=full color=0,0,255 altColor=255,255,255")
	x.hmc = cbind(x[,"chr"], x[,"pos"], x[,"pos"]+1, r.bs - r.oxbs)
	write(trackinfo, file = ofile)
	write.table(file= ofile, x.hmc, col.names=F, row.names=F, quote=F, sep="\t", append=T)


	#ratio of oxBS is estimate for mc% shown in red
	ofile = paste0(samp, "_mc.bedGraph")
	trackinfo = paste0("track type=bedGraph name=\"",ofile,"\" visibility=full altColor=255,0,0 color=255,255,255")
	x.mc = cbind(x[,"chr"], x[,"pos"], x[,"pos"]+1, -r.oxbs)
	write(trackinfo, file = ofile)
	write.table(file= ofile, x.mc, col.names=F, row.names=F, quote=F, sep="\t", append=T)

}

writeBedgraph_v2("jd", x)
writeBedgraph_v2("hpne", x)



#take a look at chr19, 15391500 - 15391600
x.19 = x[x[,"chr"] == "chr19",]

ix = (x.19[,"pos"] > 15391500) & (x.19[,"pos"] < 15391600)
x.19[ix,]


trackinfo = paste0("track type=bedGraph name=\"",bname,"\" visibility=full color=0,0,255 altColor=255,255,255")
head(x[,c("chr", "pos", "jd.bs.meth")])

#negative non-methylation red
trackinfo = paste0("track type=bedGraph name=\"",bname,"\" visibility=full altColor=255,0,0 color=255,255,255")

head(cbind(x[,c("chr", "pos")], -(x[,"jd.bs.tot"] - x[,"jd.bs.meth"])))

head(x)


