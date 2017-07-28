import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
import HTSeq as ht
import numpy.random as rd
import multiprocessing as mp

#we want to count(and eventually plot) the stat3 peaks coinciding
#with variuos region sets
#promoter regions (TSS +/- 2kb)
#genebody regions, split by exons and introns
#pancreatic regions, peak bed files specified in the ENCODE databases



#step1.  get the regions of interest
#1a.  promoters are from the annotation table file

PROM_DIST=2000
#read in the annotation file.  this was downloaded from the ucsc genome browswer
anno = np.loadtxt("/home/kpradhan/Desktop/data/hg19/anno/annotation_hg19.gtf", dtype=object)
anno[:5]
#regions lists will all be the same format
dt_sites = np.dtype([("chrom", "S4"), ("strand", "S1"), ("start", "int"), ("stop", "int")])
x = [(str(a[1]), str(a[2]), int(a[3])-PROM_DIST, int(a[3])+PROM_DIST) for a in anno]
sites_tss = np.unique(np.array(x, dtype=dt_sites))
sites_tss.shape

#1b.  can get exons/introns from the same annotaion file
#start and stop show the gene's full length
#exonStarts and exonStops will give the locations of the exons
#everything else between start and stop will be the introns
def getExonsIntrons(start, stop, exonStarts, exonStops):
    start = int(start)
    stop = int(stop)
    #split by comma
    exonStarts = map(int, exonStarts.rstrip(",").split(","))
    exonStops = map(int, exonStops.rstrip(",").split(","))
    #find all the exons
    exons = []
    for a, b in zip(exonStarts, exonStops):
        exons.append([a, b])
    #find the spots in betweent the ejdons
    introns = []
    #begin at the tss, not at the first exon
    curPos = start
    for exon in exons:
        if exon[0] - curPos > 0:
            introns.append([curPos, exon[0]])
        curPos = exon[1] 
    #remember to check the last position
    if stop - curPos > 0:
        introns.append([curPos, stop])
    return exons, introns
        


#get a list of all the exons and introns
sites_exons = []
sites_introns = []
for a in anno:
    #get the boundaries of each exon and intron
    exons, introns = getExonsIntrons(a[3], a[4], a[8], a[9])
    #add the chromosome and strand info to the list
    for e in exons:
        sites_exons.append((str(a[1]), str(a[2]), e[0], e[1]))
    for i in introns:
        sites_introns.append((str(a[1]), str(a[2]), i[0], i[1]))
#turn them into numpy record arrays
sites_exons = np.unique(np.array(sites_exons, dtype=dt_sites))
sites_introns = np.unique(np.array(sites_introns, dtype=dt_sites))
sites_exons.shape
sites_introns.shape


#1c.  load in the 6 bed files downloaded from encode
#   #H3K27ac #H3K27me3 #H3K36me3 #H3K4me1 #H3K4me3 #H3K9me3
pancFolder = "/home/kpradhan/Desktop/data/pancreas"
pancFiles = glob.glob(pancFolder+"/*.bed.gz")

pfile= pancFiles[0]
def getPeakSitesFromBed(pfile):
    peaks = np.loadtxt(pfile, dtype=object)
    #regions lists will all be the same format
    dt_sites = np.dtype([("chrom", "S5"), ("strand", "S1"), ("start", "int"), ("stop", "int")])
    x = [(str(a[1]), str(a[2]), int(a[3])-PROM_DIST, int(a[3])+PROM_DIST) for a in peaks]
    sites_tss = np.unique(np.array(x, dtype=dt_sites))


#pfile = "peaks_jd13d_stat3Input.txt"
def getPeakSitesFromMacs(pfile):
    peaks = np.loadtxt(pfile, dtype=object)
    #regions lists will all be the same format
    dt_sites = np.dtype([("chrom", "S5"), ("strand", "S1"), ("start", "int"), ("stop", "int")])
    x = []
    for a in peaks:
        #mid = (int(a[2]) + int(a[1])) / 2
        #x = [(str(a[0]), "+", mid-200, mid+200) for a in peaks]
        x.append( (str(a[0]), "+", int(a[1]), int(a[2])) )
    #only retain proper chromosomes
    x = np.unique(np.array(x, dtype=dt_sites))
    ix = np.array([a[:3] == "chr" or a == "X" or a == "Y" for a in x["chrom"]] )
    return x[ix]


#maybe only look at the peaks that pass a threshold
def getPeakMidsFromMacs(pfile, thresh=0):
    peaks = np.loadtxt(pfile, dtype=object)
    #regions lists will all be the same format
    dt_sites = np.dtype([("chrom", "S5"), ("strand", "S1"), ("site", "int")])
    x = []
    for a in peaks:
        if thresh == 0 or float(a[5]) > thresh:
            mid = (int(a[2]) + int(a[1])) / 2
            x.append((str(a[0]), "+", mid))
        #x.append( (str(a[0]), "+", int(a[1]), int(a[2])) )
    #only retain proper chromosomes
    x = np.unique(np.array(x, dtype=dt_sites))
    ix = np.array([a[:3] == "chr" or a == "X" or a == "Y" for a in x["chrom"]] )
    return x[ix]

#step2.  get the read count profiles 
#using the stat3 peaks
pfile = "peaks_jd13d_stat3Input.txt"
peaks_jd13d = getPeakSitesFromMacs(pfile)
peakmids_jd13d = getPeakMidsFromMacs(pfile)


peaks_jd13d = getPeakSitesFromMacs("peaks_jd13d_stat3Input.txt")
peaks_hpne = getPeakSitesFromMacs("peaks_hpne_stat3Input.txt")
peaks_jd13d.shape
peaks_hpne.shape


#count the number of peaks that fall in any of the regions
peaks = peaks_jd13d
sites = sites_tss
p = peaks[0]
s = sites[0]
peaks['chrom']
sites['chrom']

#peaks: X
#sites: chrX
def countOverlap(peaks, sites):
    count = 0
    #work chromosome by chrom
    for chrom in np.unique(peaks["chrom"]):
        ix1 = peaks["chrom"] == chrom
        ix2 = (sites['chrom'] == chrom) + (sites['chrom'] == "chr"+chrom)
        #does the peak region intersect any of the sites?
        A = peaks[ix1]["start"]
        B = peaks[ix1]["stop"]
        C = sites[ix2]["start"]
        D = sites[ix2]["stop"]
        for (a, b) in zip(A, B):
            count = count + np.sum((a <= D) * (C <= b))
    return count 

#find the number of peaks that lay within boundaries of
#   promoter regions
#   gene exon regions
#   gene intron regions
sites_tss.shape
sites_exons.shape
sites_introns.shape
countOverlap(peaks_jd13d, sites_tss)
countOverlap(peaks_jd13d, sites_exons)
countOverlap(peaks_jd13d, sites_introns)
countOverlap(peaks_hpne, sites_tss)
countOverlap(peaks_hpne, sites_exons)
countOverlap(peaks_hpne, sites_introns)

#for every site find the distribution of reads in
#the bam file that overlap the region centered at site
#assume the sites datastructure has a chr* name
#the bamfile might have chroms that look like "chr12" or "12"
bam = bams[0]
sites = peakmids_jd13d

class _exhausted():
    pass

def truncChr(c):
    if c.isdigit():
        return c
    elif c == "X" or c == "Y":
        return c
    else:
        return c[3:]

def getProfile(pos, bam, halfwinwidth=3000, noChr=False):
    peakProfile = np.zeros(2*halfwinwidth)
    sitechr = truncChr(pos["chrom"]) if noChr else pos["chrom"]
    window = ht.GenomicInterval( str(sitechr), pos["site"] - halfwinwidth, pos["site"] + halfwinwidth, str(pos["strand"]) )
    if (list(bam[window])):
        for almnt in bam[window]:
            if pos["strand"] == "+":
                a = almnt.iv.start - pos["site"] + halfwinwidth
                b = almnt.iv.end - pos["site"] + halfwinwidth
            if pos["strand"] == "-":
                a = pos["site"] + halfwinwidth - almnt.iv.end
                b = pos["site"] + halfwinwidth - almnt.iv.start
            peakProfile[a:b] += 1
    if (np.sum(peakProfile) > 0):
        return peakProfile
    else:
        return []

def getPeakProfiles_mt(bamfile, sites, halfwinwidth=3000, noChr = False):
    bam = ht.BAM_Reader(bamfile)
    #make sure the sites and bam files have same naming convention
    #
    #retrict sites to those that have an entry in the bam file
    bamChroms = [x["SN"] for x in bam.get_header_dict()["SQ"]]
    bamChroms = ["chr"+c if c.isdigit() else c for c in bamChroms]
    sites = sites[np.in1d(sites["chrom"], bamChroms)]
    sites.shape
    #
    peakProfs = []
    #collect the sites as genomic intervals
    pool = mp.Pool(processes=7)
    out = [pool.apply_async(getProfile, args=(pos, bam, halfwinwidth, noChr)) for pos in sites]
    peakProfs = [o.get() for o in out]
    #ignore the profiles with no reads
    peakProfs = [p for p in peakProfs if len(p) > 0]
    return np.array(peakProfs)



bamfile = "/home/kpradhan/Desktop/data/pancreas/test1/test3.bam"
def getPeakProfiles(bamfile, sites, halfwinwidth=3000, noChr = False):
    bam = ht.BAM_Reader(bamfile)
    #make sure the sites and bam files have same naming convention
    #
    #retrict sites to those that have an entry in the bam file
    bamChroms = [x["SN"] for x in bam.get_header_dict()["SQ"]]
    bamChroms = ["chr"+c if c.isdigit() else c for c in bamChroms]
    sites = sites[np.in1d(sites["chrom"], bamChroms)]
    sites.shape
    #
    peakProfs = []
    #collect the sites as genomic intervals
    for i, pos in enumerate(sites):
        peakProfile = np.zeros(2*halfwinwidth)
        if i % 1000 == 0:
            print "%d of %d" %(i, len(sites))
        sitechr = truncChr(pos["chrom"]) if noChr else pos["chrom"]
        window = ht.GenomicInterval( str(sitechr), pos["site"] - halfwinwidth, pos["site"] + halfwinwidth, str(pos["strand"]) )
        if (list(bam[window])):
            for almnt in bam[window]:
                if pos["strand"] == "+":
                    a = almnt.iv.start - pos["site"] + halfwinwidth
                    b = almnt.iv.end - pos["site"] + halfwinwidth
                if pos["strand"] == "-":
                    a = pos["site"] + halfwinwidth - almnt.iv.end
                    b = pos["site"] + halfwinwidth - almnt.iv.start
                peakProfile[a:b] += 1
        if (np.sum(peakProfile) > 0):
            peakProfs.append(peakProfile)
    return np.array(peakProfs)

def binPeakProfiles(peakProfs, K=100):
    sums = [np.sum(p) for p in peakProfs]
    ix = np.argsort(sums)
    #
    #bin every 100 rows together
    smallX = []
    for i in range(0, len(peakProfs)-1, K):
        smallX.append(np.sum(np.array([peakProfs[j] for j in ix[i:(i+K)]]),0))
    return np.array(smallX)[::-1]

peakmids_jd13d = getPeakMidsFromMacs(pfile)
sites=peakmids_jd13d
def makePeakPlots(sites, bamfile, outfile, noChr=False):
    peakProfs = getPeakProfiles(bamfile, sites, 3000, noChr)
    #
    prof = np.sum(peakProfs, 0)
    smallPeaks = binPeakProfiles(peakProfs, K=100)
    x = smallPeaks.ravel()
    vmax = np.sort(x)[np.round(len(x)*.80)]
    #
    plt.close()
    fig = plt.figure()
    fig.suptitle(outfile)
    ax = plt.subplot(121)
    xran = np.arange(len(prof))-3000
    plt.plot(xran[500:(-500)], prof[500:(-500)])
    ax = plt.subplot(122)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.imshow(smallPeaks[:, 500:(-500),],  aspect = "auto", vmax=vmax)
    #plt.colorbar(orientation="horizontal")
    plt.colorbar(orientation="vertical")
    plt.savefig(outfile)
    plt.close()


def clipFilename(ffile):
    return ffile[(ffile.rfind("/")+1):]



pfile = "peaks_jd13d_stat3Input.txt"
pancBamFolder = "/home/kpradhan/Desktop/data/pancreas/sorted/hacked"
sites_jd13d = getPeakMidsFromMacs(pfile)
sites_jd13d_q3 = getPeakMidsFromMacs(pfile, 3)
sites_jd13d_q3.shape

bams = glob.glob(pancBamFolder+"/*.bam")
for bamfile in bams:
    outfile = "jd13dStat3Peaks_q3_"+clipFilename(bamfile)+".png"
    makePeakPlots_mt(sites_jd13d_q3, bamfile, outfile)


bams = glob.glob(pancBamFolder+"/*.bam")
for bamfile in bams:
    outfile = "jd13dStat3Peaks_q3_"+clipFilename(bamfile)+".png"
    makePeakPlots(sites_jd13d, bamfile, outfile)

rsites = getRandomPeaks(sites_jd13d)
bams = glob.glob(pancBamFolder+"/*.bam")
for bamfile in bams:
    outfile = "randomPeaks_"+clipFilename(bamfile)+".png"
    makePeakPlots(rsites, bamfile, outfile)


def getRandomPeaks(sites):
    #get max sizes for each chromosome
    chroms = np.unique(sites["chrom"])
    nSites = {}
    maxSize = {}
    minSize = {}
    for chrom in chroms:
        nSites[chrom] = np.sum(sites["chrom"] == chrom)
        maxSize[chrom] = np.max(sites[sites["chrom"] == chrom]["site"])
        minSize[chrom] = np.min(sites[sites["chrom"] == chrom]["site"])
    #randomly sample peaks.
    rsites = sites.copy()
    for chrom in chroms:
        s = np.sort(rd.random_integers(minSize[chrom], maxSize[chrom], size=nSites[chrom]))
        for i, ix in enumerate(np.where(rsites["chrom"] == chrom)[0]):
            rsites[ix]["site"] = s[i]
    return rsites    
        
   
#######################################
#get the peak call from the h3k encode beds
pfile = "/media/kpradhan/DATAPART1/data/pancreas/encodePeaks/pancreas/pancreas_distal.H3K4me1.bed"
sites_d = getPeakMidsFromMacs(pfile)
pfile = "/media/kpradhan/DATAPART1/data/pancreas/encodePeaks/pancreas/pancreas_proximal.H3K4me1.bed"
sites_p = getPeakMidsFromMacs(pfile)
sites.shape
bfile = "/home/kpradhan/Desktop/data/sanchari/pancreatic/stat3_chipseq/J179/bwa/JD13D-STAT3.bwa.sorted.mkdup.bam"


pfiles = glob.glob("/media/kpradhan/DATAPART1/data/pancreas/encodePeaks/pancreas/*.bed")
bfiles = glob.glob("/home/kpradhan/Desktop/data/sanchari/pancreatic/stat3_chipseq/J179/bwa/*.bam")
for bfile in bfiles:
    for pfile in pfiles:
        sites = getPeakMidsFromMacs(pfile)
        outfile = "colocPlot_"+clipFilename(pfile)+"_"+clipFilename(bfile)+".png"
        makePeakPlots(sites, bfile, outfile, noChr=True)

#do the merged prox dist peaks
for bfile in bfiles:
    for i in range(3):
        sites = np.hstack((getPeakMidsFromMacs(pfiles[i]), getPeakMidsFromMacs(pfiles[i+3])))
        outfile = "colocPlot_"+clipFilename(pfiles[i])+clipFilename(pfiles[i+3])+"_"+clipFilename(bfile)+".png"
        makePeakPlots(sites, bfile, outfile, noChr=True)

peakProfs_p = getPeakProfiles(bfile, sites_p, 3000, noChr=True)
peakProfs_d = getPeakProfiles(bfile, sites_d, 3000, noChr=True)
prof = np.sum(peakProfs, 0)
profd = np.sum(peakProfs_d, 0)
plt.plot(prof, label="prox")
plt.plot(profd, label="distal")
plt.legend()





#################################################3
bamfile1 = "/home/kpradhan/Desktop/data/sanchari/pancreatic/stat3_chipseq/J179/bwa/jd13stat3_q10.bam"
bamfile2 = "/home/kpradhan/Desktop/data/sanchari/pancreatic/stat3_chipseq/J179/bwa/JD13D-STAT3.bwa.sorted.mkdup.bam"
bamfile2 = "/home/kpradhan/Desktop/data/sanchari/pancreatic/stat3_chipseq/J179/bwa/HPNE-Input.bwa.sorted.mkdup.bam"

sites_jd13d_p001 = getPeakMidsFromMacs(pfile, 3)
sites_jd13d = getPeakMidsFromMacs(pfile)
sites_jd13d_p001.shape
sites_jd13d.shape
peakProfs1 = getPeakProfiles(bamfile1, sites_jd13d, 3000, noChr=True)
peakProfs2 = getPeakProfiles(bamfile2, sites_jd13d, 3000, noChr=True)
prof1 = np.sum(peakProfs1, 0)
prof2 = np.sum(peakProfs2, 0)
plt.plot(prof1, label="stat3 q10")
plt.plot(prof2, label="stat3")
plt.legend()
#input has a vally
#stat3 has a mountain


bamfile = "/home/kpradhan/Desktop/data/pancreas/sorted/hacked/hacked_sorted_GSM906419_UCSD.Pancreas.Input.STL003.bam"
peakProfs = getPeakProfiles_mt(bamfile, sites_jd13d_q3, 3000)

bamfile = "/home/kpradhan/Desktop/data/pancreas/sorted/hacked/hacked_sorted_GSM906419_UCSD.Pancreas.Input.STL003.bam"
bamfile = bams[0]
peakProfs = getPeakProfiles(bamfile, sites_jd13d, 3000)
rProfs = getPeakProfiles(bamfile, rsites, 3000)

prof = np.sum(peakProfs, 0)
rprof = np.sum(rProfs, 0)
smallPeaks = binPeakProfiles(peakProfs, K=100)
rsmallPeaks = binPeakProfiles(rProfs, K=100)

plt.plot(prof)

plt.hist(smallPeaks.ravel(), 1000)
x = smallPeaks.ravel()
rx = rsmallPeaks.ravel()
vmax = np.sort(x)[np.round(len(x)*.80)]
rvmax = np.sort(rx)[np.round(len(rx)*.80)]


    plt.close()
    fig = plt.figure()
    #fig.suptitle("\n".join(map(clipFilename, [hmcfile, bamfile])))
    ax = plt.subplot(121)
    xran = np.arange(len(prof))-3000
    plt.plot(xran[500:(-500)], rprof[500:(-500)])
    #plt.plot(xran[500:(-500)], prof[500:(-500)])
    ax = plt.subplot(122)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.imshow(rsmallPeaks[:, 500:(-500),],  aspect = "auto", vmax=rvmax)
    #plt.imshow(smallPeaks[:, 500:(-500),],  aspect = "auto", vmax=vmax)
    #plt.colorbar(orientation="horizontal")
    plt.colorbar(orientation="vertical")
    #plt.savefig(outfile)
    plt.close()


np
plt.imshow(smallPeaks[:, 500:(-500),],  aspect = "auto", vmax=vmax)
plt.plot(prof[500:(-500)])
plt.plot(prof)

#is it possible to load the panc bed files as a bam or wig?
wig = ht.WiggleReader("/home/kpradhan/Desktop/data/pancreas/GSM910576_UCSD.Pancreas.H3K4me1.STL003.wig.gz")
window = ht.GenomicInterval( "chr1", 10000, 20000000, "+" )
list(wig[window])
        if (list(bam[window])):
            for almnt in bam[window]:



