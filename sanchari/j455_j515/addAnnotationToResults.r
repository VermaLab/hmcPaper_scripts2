
library("ChIPpeakAnno")
library(EnsDb.Hsapiens.v75)

library("biomaRt")
library("VariantAnnotation")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

library(GenomicFeatures)
hg19.refseq.db <- makeTxDbFromUCSC(genome="hg19", table="refGene")
refseq.genes<- genes(hg19.refseq.db)

options(stringsAsFactors=F)



getSiteCounts <- function(x){
	input = GRanges(seqnames=x[,1], ranges=IRanges(x[,2], x[,2]+1), strand="*")
	loc_hg19 <- locateVariants(input, txdb_hg19, AllVariants())
	loc_hg19 = loc_hg19[!duplicated(loc_hg19)]
	table(loc_hg19$LOCATION)
}


#read in the peak lists
atac.hpne.only = read.table("macs_Cjd_sorted_mdup_merged_Thpne_sorted_mdup/Cjd_sorted_mdup_merged_Thpne_sorted_mdup_peaks.xls", header=T)
atac.jd.only = read.table("macs_Chpne_sorted_mdup_Tjd_sorted_mdup_merged/Chpne_sorted_mdup_Tjd_sorted_mdup_merged_peaks.xls", header=T)
#atac.hpne.only = read.table("macs_Cjd_sorted_mdup_merged_Thpne_sorted_mdup/Cjd_sorted_mdup_merged_Thpne_sorted_mdup_summits.bed")
#atac.jd.only = read.table("macs_Chpne_sorted_mdup_Tjd_sorted_mdup_merged/Chpne_sorted_mdup_Tjd_sorted_mdup_merged_summits.bed")
head(atac.hpne.only)
head(atac.jd.only)

#this is how it's done with oxbs
sCounts.hpne = getSiteCounts(atac.hpne.only)
sCounts.jd = getSiteCounts(atac.jd.only)


sCounts.hpne 
sCounts.jd 

sCounts.hpne/sum(sCounts.hpne)
sCounts.jd /sum(sCounts.jd)





chunk <- function(x, n){
	split(x, ceiling(seq_along(x)/n))
}

getGenesFromPeaks_chunks <- function(peakValue_chunks){
	results = list()
	#process in chunks
	for (i in 1:length(peakValue_chunks)){
		chunk = peakValue_chunks[[i]]
		print (chunk[1])
		res = NA
		try({
			G = getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region"), values = chunk, mart=grch37)
			res = unique(G[,'hgnc_symbol'])
		}, silent = T)
		results[[i]] = res
	}

	results
}


getGenesFromPeaks <- function(peaks, chunkSize=1000){
	values = paste0(gsub(peaks[,1], pat="chr(.*)", rep="\\1", perl=T),":", peaks[,2], ":", peaks[,3])
	peakValue_chunks = chunk(values, chunkSize)

	results = getGenesFromPeaks_chunks(peakValue_chunks)

	unique(unlist(results))
	#or maybe do a while loop
	#if there were any errors, try one more time
	#check results for NA entries
	#ix = which(sapply(results, function(x){any(is.na(x))}))
	#if (length(ix) > 0){
	#	results2 = getGenesFromPeaks_chunks(peakValue_chunks[ix])
	#}

}

filterPeaks <- function(peaks, grange){
	#intersect the peaks with the TSS locations
	#the chromosomal locations of the peaks
	#granges requires format of chr1:pos1-pos2
	x = GRanges(paste0(gsub(peaks[,1], pat="chr(.*)", rep="\\1", perl=T),":", peaks[,2], "-", peaks[,2]+1))
	#which peaks overlap with Tss's
	res = findOverlaps(x, grange)
	ix = unique(queryHits(res))
	peaks[ix,]
}

#get a list of TSS's +/- bp for all all genes
all.genes <- unique( getBM(attributes = "hgnc_symbol", values = "*", mart = grch37) )
bp = 2000
tss = getBM(attributes=c("chromosome_name", "transcript_start", "transcript_end"),
      filters="hgnc_symbol", values=all.genes$hgnc_symbol, mart=grch37)
tss.2k =  GRanges(paste0(tss$chromosome_name,":", tss$transcript_start-bp, "-", tss$transcript_start+bp))



#filter away all the sites that are more than 2kb from a TSS
atac.jd.tss2k = filterPeaks(atac.jd.only, tss.2k)
atac.hpne.tss2k = filterPeaks(atac.hpne.only, tss.2k)
dim(atac.jd.tss2k)
dim(atac.hpne.tss1k)



appendClosestGeneInfo <- function(x, pos = 2){
	res = lapply(1:nrow(x), function(i){
		if (i %% 1000 == 0){
			print(i)
		}
		genes = findClosestGene(x[i,1],x[i, pos],"hg19")
		#if there's more than one, just take the first
		data.frame(genes[1,])
	})

	x.res = do.call(rbind, res)
	cbind(x.res, x)
}

atac.jd.tss2k.genes = appendClosestGeneInfo(atac.jd.tss2k, 5)
atac.hpne.tss2k.genes = appendClosestGeneInfo(atac.hpne.tss2k, 5)

head(atac.jd.tss2k.genes)
head(atac.hpne.tss2k.genes)

#filter away the peaks further than 2k from a gene tss
filterOnTss <- function(x, D){
	x[abs(x$abs_summit - x$txStart) <= D,]
}


atac.jd.tss2k.genes[abs(atac.jd.tss2k.genes$Distance) <= 2000,]


x.jd  = filterOnTss(atac.jd.tss2k.genes, 2000)
x.hpne = filterOnTss(atac.hpne.tss2k.genes, 2000)
write.table(x = x.jd, file="geneAnnoSets/atac_jd_tss2k_genes.txt", row.names=F, sep="\t")
write.table(x = x.hpne, file="geneAnnoSets/atac_hpne_tss2k_genes.txt", row.names=F, sep="\t")


#redo the region type counts
#on just the TSS close peaks
sCounts.hpne = getSiteCounts(x.jd)
sCounts.jd = getSiteCounts(x.hpne)

sCounts.hpne 
sCounts.jd 

sCounts.hpne/sum(sCounts.hpne)
sCounts.jd /sum(sCounts.jd)



#######################################


x = atac.jd.tss2k


res = lapply(1:nrow(x), function(i){
	if (i %% 1000 == 0){
		print(i)
	}
	genes = findClosestGene(x[i,1],x[i, 2],"hg19")
	#if there's more than one, just take the first
	data.frame(genes[1,])
})

x.res = do.call(rbind, res)
head(x.res)



#load in all the gene info from hg19
knownGenes <- genes(EnsDb.Hsapiens.v75)
rd <- RangedData(IRanges(start = x[,2],
	end = x[,3]), space = x[,1])
ranges.peaks = toGRanges(rd, format="RangedData")

#annotate ranges with the knownGenes database
anno <- annotatePeakInBatch(ranges.peaks, AnnotationData=knownGenes)


#make sure to only take 1 result
x = anno[!(duplicated(anno$peak))]
x$feature






tssgenes.jd = getGenesFromPeaks(atac.jd.tss2k)
tssgenes.hpne = getGenesFromPeaks(atac.hpne.tss2k)








lis

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
head(atac.hpne.only)
pos = paste0(substr(atac.hpne.only[,1], 4, 99), ":", atac.hpne.only[,2], ":", atac.hpne.only[,3])
#how to get gene id info
G = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype"), filters = "chromosomal_region", values = pos[1:1000], mart=grch37)
