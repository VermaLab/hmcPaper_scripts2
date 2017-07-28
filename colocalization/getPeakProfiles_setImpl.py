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



#$reload(pp)


outFolder = "test4"
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


#create a genomic array of sets to hold the site info
gas = ht.GenomicArrayOfSets("auto", stranded=True)
?ht.GenomicInterval()

def site2GI(site):
    return ht.GenomicInterval(""+site[0], site[2]-W, site[2]+W, ""+site[1])

#fill the arary of sets
for i, site in enumerate(sites[:10]):
    iv = site2GI(site)
    gas[iv] += str(i)

for g in gas.steps():
    print g


gas['1']

for iv, val in gas[ iv1 ].steps():
    print iv, val   

site = sites[1]
iv1 = ht.GenomicInterval(""+site[0], site[2]-W, site[2]+W, ""+site[1])

#get regions spanned by site 1
site-W, 


#read the bam file read by read filling in genomic array of sets


#look at coverages of sites















