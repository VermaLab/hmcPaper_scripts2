#from multiprocessing import Pool
from functools import partial
from itertools import chain
import time
import peakProfiles as pp
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import HTSeq as ht
import random



#reload(pp)



outFolder = "test5"
#siteFile  = "hpneHmc_strandedPeakSummits.txt"
siteFile  = "hpne_mc.txt"
siteFile  = "hpne_stat3.txt"
siteFile = "hpne_hmc_p05.bed"
#bamfile1  = "/home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"
#bamfile1 = "/home/kpradhan/Desktop/hpc_home/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"
bamfile1 = "/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/hpne_sorted_mdup.bam"
bamfile1 = "/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455/hpne_sorted_mdup_chr10.bam"
W         = 3000
#K=200000  #only look at this number of sites due to memory constraints
K=400000  #only look at this number of sites due to memory constraints

#command line arguments
outFolder = sys.argv[1]
siteFile  = sys.argv[2]
bamfile1  = sys.argv[3]
W         = int(sys.argv[4])

print ("outfolder:  %s" % outFolder)
print ("siteFile:  %s" % siteFile)
print ("bamfile:  %s" % bamfile1)
print ("W:  %d" % W)



#load the sites
sites = pp.getPeakSitesFromStrandedBed(siteFile)


#turn any sites with * strands into .
for s in sites:
    if s['strand'] != "+" and s['strand'] != "-":
        s['strand'] = "."


#take out any sites that are lower than the halfwidth
ix = sites['site'] > W
sites = sites[ix]

#first thing...check that the sites match the notation of the bam file
sites, nMatch = pp.checkSiteNotation(sites, bamfile1)
len(sites)
print("number of sites matching bam chroms: %d" % nMatch)


#get rid of the sites that have no chromosome entry in the bam file

#save results in a folder
if not os.path.exists(outFolder):
    os.makedirs(outFolder)

#save some batch information
with open(os.path.join(outFolder, "batch_info.txt"), "w") as myfile:
    myfile.write("outfolder:  %s\n" % outFolder)
    myfile.write("siteFile:  %s\n" % siteFile)
    myfile.write("bamfile:  %s\n" % bamfile1)
    myfile.write("W:  %d\n" % W)
    myfile.write("#matching sites:  %d\n" % nMatch)


#if all the sites are facing the same way, used * strand coverage
stranded = True
if len(np.unique(sites['strand'])) == 1:
    stranded = False
    sites['strand'] = "."


#load the entire bam into a genomic array
bam = ht.BAM_Reader(bamfile1)
#maybe, if all the site strands are *, don't strand coverage
coverage = ht.GenomicArray( "auto", stranded=stranded, typecode="i" )
#coverage = ht.GenomicArray( "auto", stranded=True, typecode="i" )
for i, almnt in enumerate(bam):
    if i % 100000 == 0:
        print (i )
    #check for quality and pcr dups!
    if almnt.aligned and almnt.aQual > 0 and not almnt.pcr_or_optical_duplicate:
        coverage[ almnt.iv ] += 1


#reload(pp)
#todo
##just save the summed Profile
peakProf = pp.getProfileFromCoverage(sites, coverage, halfwinwidth=3000, stranded=stranded)

np.save(os.path.join(outFolder, "peakSums.npy"), peakProf)



