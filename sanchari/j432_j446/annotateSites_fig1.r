library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) # for annotation
library(org.Hs.eg.db)   
library("biomaRt")
library("stringr")

options(stringsAsFactors=F)

#load up the database utlility, find the chrom pos and ids for the significant diff expressed locatinos 
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
#grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")






#https://adairama.wordpress.com/2013/02/15/functionally-annotate-snps-and-indels-in-bioconductor/
annotateList <- function(sites){

	#get the chromosomal location of the sites
	input = sites[, 1:2]
	colnames(input) <- c("chr", "pos")
	input$pos       <- as.numeric(as.character(input$pos))
	 
	#construct a range for each site.
	target <- with(input, GRanges( 
		seqnames = Rle(chr),
		ranges   = IRanges(
			pos, 
		    end=pos, 
		    names=paste0("rand", 1:nrow(input)) 
		),
		strand   = Rle(strand("*"))
	))

	#find variants for each site
	loc <- locateVariants(target, TxDb.Hsapiens.UCSC.hg19.knownGene, AllVariants())
	names(loc) <- NULL
	out <- as.data.frame(loc)
	out$names <- names(target)[ out$QUERYID ]
	out <- out[ , c("names", "seqnames", "start", "end", "LOCATION", "GENEID")]
	out <- unique(out)

	#convert gene id to gene symbol
	Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
	id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
	names(id2Symbol) <- unlist(Symbol2id)
	 
	#make sure all entries have been found
	x <- unique( with(out, c(levels(GENEID))) )
	table( x %in% names(id2Symbol) ) # good, all found
	 
	#assign the gene symbol
	out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
	out
}


makePie <- function(x, descr){
	mytable <- table(x[,5])
	#lbls <- paste(names(mytable), "\n", mytable, sep="")
	lbls <- paste(mytable, " (", format(round(100*mytable/sum(mytable), 4), nsmall=4), "%)", sep="")
	cols = rainbow(length(lbls))
	pie(mytable, labels = lbls, col=cols, main=descr) 
	legend("topright", names(mytable), col=cols, pch=15)
	mytable
}



#load in the sites
fol = "/home/kpradhan/Desktop/projects/sanchari/homer_20160708"
atac.hpne = read.table(paste0(fol, "/atac_hpne.bed"), header=F)
atac.jd = read.table(paste0(fol, "/atac_jd.bed"), header=F)
atac.hpne.only = read.table(paste0(fol, "/atac_hpne_only.bed"), header=F)
atac.jd.only = read.table(paste0(fol, "/atac_jd_only.bed"), header=F)

hmc.hpne.only = read.table(paste0(fol, "/hmc_hpne_only.bed"), header=F)
hmc.jd.only = read.table(paste0(fol, "/hmc_jd_only.bed"), header=F)

mc.hpne.only = read.table(paste0(fol, "/mc_hpne_only.bed"), header=F)
mc.jd.only = read.table(paste0(fol, "/mc_jd_only.bed"), header=F)


atac.hpne.anno = annotateList(atac.hpne)
atac.jd.anno = annotateList(atac.jd)

atac.hpne.only.anno = annotateList(atac.hpne.only)
atac.jd.only.anno = annotateList(atac.jd.only)
hmc.hpne.only.anno = annotateList(hmc.hpne.only)
hmc.jd.only.anno = annotateList(hmc.jd.only)
mc.hpne.only.anno = annotateList(mc.hpne.only)
mc.jd.only.anno = annotateList(mc.jd.only)

png("pie_atac_hpne.png")
makePie(atac.hpne.anno, "atac_hpne")
dev.off()
png("pie_atac_jd.png")
makePie(atac.jd.anno, "atac_jd")
dev.off()

png("pie_atac_hpne_only.png")
makePie(atac.hpne.only.anno, "atac_hpne_only")
dev.off()
png("pie_atac_jd_only.png")
makePie(atac.jd.only.anno, "atac_jd_only")
dev.off()

png("pie_hmc_hpne_only.png")
makePie(hmc.hpne.only.anno, "hmc_hpne_only")
dev.off()
png("pie_hmc_jd_only.png")
makePie(hmc.jd.only.anno, "hmc_jd_only")
dev.off()

png("pie_mc_hpne_only.png")
makePie(mc.hpne.only.anno, "mc_hpne_only")
dev.off()
png("pie_mc_jd_only.png")
makePie(mc.jd.only.anno, "mc_jd_only")
dev.off()














#look at a random distribution of sites taken from min/max range of chromosome
sites.rand = sites.jd
max(sites.rand[,2])
sites.rand[,2] = runif(
	nrow(sites.rand),
	min(sites.rand[,2]), 
	max(sites.rand[,2]) - min(sites.rand[,2])
)




#break this up...use just 10% of the sites
#to speed things up
#annoation the mcSites
#break up long list into multiple queries
#takes a long time
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))
P = 0.1
ix.small = sort(sample(1:nrow(sites.hpne), nrow(sites.hpne)*P))
smallsites.hpne = sites.hpne[ix.small,]
ix = chunk(1:nrow(smallsites.hpne), 10)
anno.hpne = lapply(ix, function(i){
	print (i[1])
	annotateList(smallsites.hpne[i,])
})
x.hpne=do.call(rbind, anno.hpne)

ix.small = sort(sample(1:nrow(sites.jd), nrow(sites.jd)*P))
smallsites.jd = sites.jd[ix.small,]
ix = chunk(1:nrow(smallsites.jd), 10)
anno.jd = lapply(ix, function(i){
	print (i[1])
	annotateList(smallsites.jd[i,])
})
x.jd=do.call(rbind, anno.jd)

save.image("mc_anno_v1.RData")


png("pie_mc_Hpne_all.png")
makePie(x.hpne, "HPNE mc sites")
dev.off()
png("pie_mc_JD_all.png")
makePie(x.jd, "JD mc sites")
dev.off()

#how many genes are affected?
head(x.jd)
ugenes.jd = unique(x.jd$GENESYMBOL)
length(ugenes.jd)
ugenes.hpne = unique(x.hpne$GENESYMBOL)
length(ugenes.hpne)


#genes with hmc in one sample, but not the other
length(ugenes.hpne[!(ugenes.hpne %in% ugenes.jd)])
length(ugenes.jd[!(ugenes.jd %in% ugenes.hpne)])





anno.hpne = annotateList(sites.hpne)
anno.jd = annotateList(sites.jd)
anno.rand = annotateList(sites.rand)
anno.hpne.sig = annotateList(sites.hpne.sig)
anno.jd.sig = annotateList(sites.jd.sig)
dim(sites.hpne)
dim(sites.jd)
dim(sites.hpne.sig)
dim(sites.jd.sig)


#briefly look at hmc events at brd4
#save in bed format, view in igv
ix.1 = which(anno.hpne$GENESYMBOL == "BRD4")
ix.2 = which(anno.jd$GENESYMBOL == "BRD4")
length(ix.1)
length(ix.2)

anno.hpne[ix.1,]
x = anno.hpne[ix.1,2:4]
x = cbind(x, 1, 1)
write.table(file="hmc_brd4_hpne.bed", x=x, col.names=F, row.names=F, quote=F)
x = anno.jd[ix.2,2:4]
x = cbind(x, 1, 1)
write.table(file="hmc_brd4_jd.bed", x=x, col.names=F, row.names=F, quote=F)
#add a header to each bed file
# track name=test1 description="hmc" visibility=2


png("atacseq_hpne.png")
mytable <- table(panno.hpne[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="ATACseq: HPNE") 
dev.off()

png("atacseq_jd.png")
mytable <- table(panno.jd[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="ATACseq: JD") 
dev.off()


#par(mfcol=c(1, 2))
png("hmc_chr19_hpne.png")
mytable <- table(anno.hpne[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="chr19 HMC: HPNE") 
dev.off()

png("hmc_chr19_jd.png")
mytable <- table(anno.jd[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="chr19 HMC: JD") 
dev.off()

png("hmc_chr19_hpneSig.png")
mytable <- table(anno.hpne.sig[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="chr19 HMC: HPNEsig") 
dev.off()

png("hmc_chr19_jdSig.png")
mytable <- table(anno.jd.sig[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="chr19 HMC: JDsig") 
dev.off()

png("hmc_chr19_random.png")
mytable <- table(anno.rand[,5])
lbls <- paste(names(mytable), "\n", mytable, sep="")
pie(mytable, labels = lbls,
   main="chr19 HMC: random") 
dev.off()

#doesn't look like there's much differenc in percentage
#on exonic/intron/promoter regions between hpne, jd
#random sites have slighter more intergenic regions


#might have to filter out some of the hmc sites based
#on direction of the gene
#wehther the reference base is C->T or G->A


#maybe look at the mutation impact factors?
#with snpEff?
#or with VEP

#Lets look at the hmc sites at brd4
#prepare teh vep file

sites = anno.hpne[ix.1,2:4]


anno.hpne[ix.1,]
x = anno.hpne[ix.1,2:4]



x = cbind(x, 1, 1)
write.table(file="hmc_brd4_hpne.bed", x=x, col.names=F, row.names=F, quote=F)
x = anno.jd[ix.2,2:4]


dat = GRanges(seqnames=rep('chr17',3),ranges=IRanges(start=c(45229228,45229234,45234706),width=1))
varallele = DNAStringSet(c('C','C','C'))
db = TxDb.Hsapiens.UCSC.hg19.knownGene

results = predictCoding(dat,db,Hsapiens,varallele)


sites = sites.hpne
input = sites[, 1:2]
colnames(input) <- c("chr", "pos")
input$pos       <- as.numeric(as.character(input$pos))
input = input[1:3000,] 

target <- with(input, GRanges( 
	seqnames = Rle(chr),
	ranges   = IRanges(
		pos, 
		width=1, 
		names=paste0("rand", 1:nrow(input)) 
	),
	strand   = Rle(strand("*"))
))
length(target)
varallele = DNAStringSet(rep('C',length(target)))
results = predictCoding(target,db,Hsapiens,varallele)






#############################################
"chromosome_name", "start", "end"
x = list(chromosome_name="19", start="15236836", end="15332543")

x = head(sites.hpne[,c(1,2,2)], 300)
x = list(chromosome_name=gsub(x[,1], pat="chr(.*)", rep="\\1", perl=T), start=x[,2], end=x[,3])
x = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "go_biological_process_linkage_type"), filters = c("chromosome_name", "start", "end"), values = x, mart=grch37)

