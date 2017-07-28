library("ChIPpeakAnno")
#source("https://bioconductor.org/biocLite.R")
#biocLite("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)

options(stringsAsFactors=F)

#load in all the gene info from hg19
knownGenes <- genes(EnsDb.Hsapiens.v75)




#peaksFile must be the peak summit file from macs2
peaksFile = "hpne_stat3_summits.bed"
addStrandToPeaks <- function(peaksFile, newFile){
	#read in summits and convert to ranges
	peaks = read.table(peaksFile)
	#if the 3rd column is not the peak ends
	if (is.na(as.numeric(peaks[1,3]))){
		#use the peaks starts + 1
		ends = peaks[,2] + 1
	}else{
		ends = peaks[,3]
	}
	rd <- RangedData(IRanges(start = peaks[,2],
		end = ends), space = peaks[,1])
	ranges.peaks = toGRanges(rd, format="RangedData")

	#annotate ranges with the knownGenes database
	anno <- annotatePeakInBatch(ranges.peaks, AnnotationData=knownGenes)

	#make sure to only take 1 result
	x = anno[!(duplicated(anno$peak))]
	#append the strand info to last column
	y = cbind(peaks, strand=as.character(x$feature_strand))
	y$strand = as.character(x$feature_strand)
	y$strand[is.na(y$strand)] = "*"
	#write the new stranded peak file
	write.table(newFile, x=y, col.names=F, row.names=F, quote=F, sep="\t")
}


#stat3
addStrandToPeaks("test1_jdmc_top100k.txt", "test1_jdmc_top100k_stranded.txt")
addStrandToPeaks("test1_jdhmc.txt", "test1_jdhmc_stranded.txt")

addStrandToPeaks("hpne_stat3_summits.bed", "hpne_stat3_summits_stranded.txt")
addStrandToPeaks("jd13d_stat3_summits.bed", "jd_stat3_summits_stranded.txt")

peaksFile = "peaks.atac.jd.gain.1k.txt"
addStrandToPeaks("peaks.atac.jd.gain.1k.txt", "peaks.atac.jd.gain.1k_stranded.txt")
addStrandToPeaks("peaks.atac.jd.loss.1k.txt", "peaks.atac.jd.loss.1k_stranded.txt")


addStrandToPeaks("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_Chpne_sorted_mdup_Tjd_sorted_mdup_merged/Chpne_sorted_mdup_Tjd_sorted_mdup_merged_summits.bed", "jd-hpne_strandedPeakSummits.txt")
addStrandToPeaks("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_Cjd_sorted_mdup_merged_Thpne_sorted_mdup/Cjd_sorted_mdup_merged_Thpne_sorted_mdup_summits.bed", "hpne-jd_strandedPeakSummits.txt")


addStrandToPeaks("p3.txt", "p3_stranded.txt")

addStrandToPeaks("hmc_sig0.001_commonSites1_hpne.bed", "hpneHmc_strandedPeakSummits_0.001.txt")
addStrandToPeaks("hmc_sig0.001_commonSites1_jd.bed", "jdHmc_strandedPeakSummits_0.001.txt")

addStrandToPeaks("hmc_sig0.01_commonSites1_hpne.bed", "hpneHmc_strandedPeakSummits.txt")
addStrandToPeaks("hmc_sig0.01_commonSites1_jd.bed", "jdHmc_strandedPeakSummits.txt")

addStrandToPeaks("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/macs_hpne_sorted_mdup/hpne_sorted_mdup_summits.bed", "hpne_strandedPeakSummits.txt")
addStrandToPeaks("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_jd_sorted_mdup_merged/jd_sorted_mdup_merged_summits.bed", "jd_strandedPeakSummits.txt")


addStrandToPeaks("/home/kpradhan/mnt/hpc_home/projects/colocalization/jd_hmc.txt", "jd_strandedHmc.txt")
addStrandToPeaks("/home/kpradhan/mnt/hpc_home/projects/colocalization/hpne_hmc.txt", "hpne_strandedHmc.txt")



#######################################################33
#scratch









rd <- RangedData(IRanges(start = peaks.hpne[,2],
	end = peaks.hpne[,3]), space = peaks.hpne[,1])
ranges.hpne = toGRanges(rd, format="RangedData")

anno.hpne <- annotatePeakInBatch(ranges.hpne, AnnotationData=knownGenes)

x = anno.hpne[!(duplicated(anno.hpne$peak))]
#write the peaks with strand info
x$feature_strand
slot(strand(x), "values")
str(ranges(x))
y = cbind(peaks.hpne, strand=as.character(x$feature_strand))
y$strand = as.character(y$strand)
class(y$strand)
sum(is.na(y$strand))
y$strand[is.na(y$strand)] = "*"

write.table("peaks_test1.bed", x=y, col.names=F, row.names=F, quote=F, sep="\t")

head(y)

start(ranges(x))
end(ranges(x))
head(x)
str(x)
names(x)
x$ranges
x["ranges"]
ranges(x)
dim(peaks.hpne)



hmc_sig0.01_commonSites1_hpne.bed

hmc.hpne = read.table("hmc_sig0.01_commonSites1_hpne.bed")

peaks.hpne = read.table("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/macs_hpne_sorted_mdup/hpne_sorted_mdup_summits.bed")
peaks.hpne = read.table("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/macs_hpne_sorted_mdup/hpne_sorted_mdup_peaks.xls", header=T)
head(peaks.hpne)

peaks.jd = read.table("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_jd_sorted_mdup_merged/jd_sorted_mdup_merged_summits.bed")

peaks.jdMinushpne = read.table("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_Chpne_sorted_mdup_Tjd_sorted_mdup_merged/Chpne_sorted_mdup_Tjd_sorted_mdup_merged_summits.bed")
peaks.hpneMinusjd = read.table("/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/macs_Cjd_sorted_mdup_merged_Thpne_sorted_mdup/Cjd_sorted_mdup_merged_Thpne_sorted_mdup_summits.bed")

