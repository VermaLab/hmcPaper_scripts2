import numpy as np

import matplotlib.pyplot as plt
import glob
import os
import sys
import HTSeq as ht
import numpy.random as rd


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


#the bed file should have 3 columns
#chrom
#position start
#position end
#strand + or -
def getPeakSitesFromStrandedBed(pfile):
    peaks = np.loadtxt(pfile, dtype=object)
    #regions lists will all be the same format
    dt_sites = np.dtype([("chrom", "S5"), ("strand", "S1"), ("site", "int")])
    #x = [(str(a[0]), str(a[3]), int(a[2])) for a in peaks]
    x = [(str(a[0]), str(a[2]), int(a[1])) for a in peaks]
    return np.unique(np.array(x, dtype=dt_sites))



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


def truncChr(c):
    if c.isdigit():
        return c
    elif c == "X" or c == "Y":
        return c
    else:
        return c[3:]


#how many sites have thee same notation as the bam file.
def checkSiteNotation(sites, bamfile):
    bam = ht.BAM_Reader(bamfile)
    #make sure the sites and bam files have same naming convention
    #
    #retrict sites to those that have an entry in the bam file
    #probably has bug dealing with X and Y.
    bamChroms = [x["SN"] for x in bam.get_header_dict()["SQ"]]
    sites = sites[np.in1d(sites["chrom"], bamChroms)]
    return sites, np.sum(np.in1d(sites["chrom"], bamChroms))


def getProfileFromCoverage(sites, coverage, halfwinwidth=3000, stranded=False):
    nSites = 0
    peakSums = np.zeros(2*halfwinwidth)
    #collect the sites as genomic intervals
    for i, pos in enumerate(sites):
        #print "%d of %d" %(i, len(sites))
        peakProfile = np.zeros(2*halfwinwidth)
        if i % 1000 == 0:
            print "%d of %d" %(i, len(sites))
        leftSide = pos["site"] - halfwinwidth
        if leftSide < 1:
            leftSide = 1
        rightSide = pos["site"] + halfwinwidth
        if (stranded):
            window = ht.GenomicInterval( str(pos["chrom"]), leftSide, rightSide, str(pos["strand"]) )
        else:
            window = ht.GenomicInterval( str(pos["chrom"]), leftSide, rightSide, ".")
        #print (pos)
        #try this out 1-17-17
        wincvg = np.array(list(coverage[window]))
        #if the iterator has nothing, it'll throw an error
        #wincvg = np.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
        if (stranded):
            if pos["strand"] == "+":
                peakProfile += wincvg
            elif pos["strand"] == "-":
                peakProfile += wincvg[::-1]
        else:
            peakProfile += wincvg
        if (np.sum(peakProfile) > 0):
            peakSums += peakProfile
    return peakSums


def getPeakProfilesFromCoverage(sites, coverage, halfwinwidth=3000, stranded=False):
    peakProfs = []
    #collect the sites as genomic intervals
    for i, pos in enumerate(sites):
        #print "%d of %d" %(i, len(sites))
        peakProfile = np.zeros(2*halfwinwidth)
        if i % 1000 == 0:
            print "%d of %d" %(i, len(sites))
        leftSide = pos["site"] - halfwinwidth
        if leftSide < 1:
            leftSide = 1
        rightSide = pos["site"] + halfwinwidth
        if (stranded):
            window = ht.GenomicInterval( str(pos["chrom"]), leftSide, rightSide, str(pos["strand"]) )
        else:
            window = ht.GenomicInterval( str(pos["chrom"]), leftSide, rightSide, ".")
        #print (pos)
        #try this out 1-17-17
        wincvg = np.array(list(coverage[window]))
        #if the iterator has nothing, it'll throw an error
        #wincvg = np.fromiter( coverage[window], dtype='i', count=2*halfwinwidth )
        if (stranded):
            if pos["strand"] == "+":
                peakProfile += wincvg
            elif pos["strand"] == "-":
                peakProfile += wincvg[::-1]
        else:
            peakProfile += wincvg
        if (np.sum(peakProfile) > 0):
            peakProfs.append(peakProfile)
    return np.array(peakProfs)



#need a way to determine if the bam file has notation as "chr10" or "10"

bamfile = "/home/kpradhan/Desktop/data/pancreas/test1/test3.bam"
def getPeakProfiles(sites, bamfile,halfwinwidth=3000):
    bam = ht.BAM_Reader(bamfile)
    #make sure the sites and bam files have same naming convention
    #
    #retrict sites to those that have an entry in the bam file
    #probably has bug dealing with X and Y.
    bamChroms = [x["SN"] for x in bam.get_header_dict()["SQ"]]
    bamChroms = ["chr"+c if c.isdigit() else c for c in bamChroms]
    sites = sites[np.in1d(sites["chrom"], bamChroms)]
    sites.shape
    #
    peakProfs = []
    #collect the sites as genomic intervals
    for i, pos in enumerate(sites):
        print "%d of %d" %(i, len(sites))
        peakProfile = np.zeros(2*halfwinwidth)
        if i % 1000 == 0:
            print "%d of %d" %(i, len(sites))
        #don't change the site notation here
        #sitechr = truncChr(pos["chrom"]) if hasChrPrefix else pos["chrom"]
        sitechr = pos["chrom"]
        window = ht.GenomicInterval( str(sitechr), pos["site"] - halfwinwidth, pos["site"] + halfwinwidth, str(pos["strand"]) )
        #if (list(bam[window])):
        if next(bam[window], None) is not None:
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


#prof is the aveage profile
#smallPeaks is the binned peak image
#cuts off 500 from either side
def savePeakPlots(outfile, prof, smallPeaks):
    #half width of the window
    W = len(prof)/2
    x = smallPeaks.ravel()
    vmax = np.sort(x)[np.round(len(x)*.80)]
    #
    plt.close()
    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    fig.suptitle(outfile)
    ax = plt.subplot(121)
    xran = np.arange(len(prof))-W
    plt.plot(xran[500:(-500)], prof[500:(-500)])
    ax = plt.subplot(122)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.imshow(smallPeaks[:, 500:(-500),],  aspect = "auto", vmax=vmax)
    #plt.colorbar(orientation="horizontal")
    plt.colorbar(orientation="vertical")
    plt.savefig(outfile, dpi=100)
    plt.close()



def savePeakPlots_part(outfile, prof, smallPeaks, throwPortion=0.20):
    rowStart = int(smallPeaks.shape[0]*throwPortion )
    rowEnd = smallPeaks.shape[0]
    #
    #half width of the window
    W = len(prof)/2
    x = smallPeaks.ravel()
    vmax = np.sort(x)[np.round(len(x)*.80)]
    #
    plt.close()
    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    fig.suptitle(outfile)
    ax = plt.subplot(121)
    xran = np.arange(len(prof))-W
    #
    test1 = np.sum(smallPeaks[rowStart:rowEnd, 500:(-500),],0)
    plt.plot(xran[500:(-500)], test1)
    np.sum(smallPeaks[1000:4000, 500:(-500),],0)
    #
    #plt.plot(xran[500:(-500)], prof[500:(-500)])
    ax = plt.subplot(122)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    #
    x = smallPeaks[rowStart:rowEnd, 500:(-500),].ravel()
    vmax = np.sort(x)[np.round(len(x)*.80)]
    plt.imshow(smallPeaks[rowStart:rowEnd, 500:(-500),],  aspect = "auto", vmax=vmax)
    #plt.colorbar(orientation="horizontal")
    plt.colorbar(orientation="vertical")
    plt.savefig(outfile, dpi=100)
    plt.close()


#test


#peakmids_jd13d = getPeakMidsFromMacs(pfile)
#sites=peakmids_jd13d
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


