#from multiprocessing import Pool
from functools import partial
from itertools import chain
import time
import peakProfiles as pp
import numpy as np
import matplotlib.pyplot as plt
import os
import sys



#$reload(pp)


outFolder = "test3"
#siteFile  = "hpneHmc_strandedPeakSummits.txt"
siteFile  = "hpne_atac.txt"
#bamfile1  = "/home/kpradha1/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"
bamfile1 = "/home/kpradhan/Desktop/hpc_home/projects/sanchari/rnaseq/HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam"
W         = 3000

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
#first thing...check that the sites match the notation of the bam file
nMatch = pp.checkSiteNotation(sites, bamfile1)
print("number of sites matching bam chroms: %d" % nMatch)

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


#partial function that allows multiple arguments to be passed to map
getPPhelper = partial(pp.getPeakProfiles, bamfile = bamfile1, halfwinwidth=W)

#pp.getPeakProfiles(sites[:100],bamfile1, W, hasChrPrefix=False)
#getPPhelper(sites[:100])


#multiprocessing
# Make the Pool of workers
#pool = Pool(processes=4)


b = time.time()
print (b)
#results = pool.map(getPPhelper, pp.chunks(sites, 1000))
#results = pool.map(getPPhelper, pp.chunks(sites[:100], 10))
results = pp.getPeakProfiles(sites[:100,  bamfile = bamfile1, halfwinwidth=W)
#results = getPPhelper(sites[:100])
peakProfs = np.concatenate(results)
a = time.time()
print (a - b)

#close the pool and wait for the work to finish (neccesary?)
#pool.close() 
#pool.join() 


#or mabye use the mean
prof = np.sum(peakProfs, 0)
#the binned image
smallPeaks = pp.binPeakProfiles(peakProfs, K=100)

np.save(os.path.join(outFolder, "peakProfiles.npy"), peakProfs)
np.save(os.path.join(outFolder, "summedPP.npy"), prof)
np.save(os.path.join(outFolder, "binnedPP.npy"), smallPeaks)





