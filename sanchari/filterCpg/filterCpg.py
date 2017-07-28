#!/usr/bin/python
import numpy as np
from itertools import islice
import os
import sys

faFol = "/home/kpradhan/Desktop/hpc_home/projects/hg19/chroms"
faFol = "/home/kpradha1/projects/hg19/chroms"
posFile = "/home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/commonSites/common4.txt"
outFile = "test1.txt"

def loadFasta(faFile):
    with open(faFile) as fp:
        fa = ""
        #load the first line which should start with >
        line = fp.readline()
        #print (line)
        fa = "".join([line.rstrip("\n") for line in fp.readlines()])
        return fa 




#open the output file
#out = open(outFile, "w")

curChrom = ""
#for line in islice(open(posFile), 1000):
#process chromosomal positions line by line
#for line in open(posFile):
for line in sys.stdin:
    line = line.rstrip("\n")
    #print line
    #first two fields are the chromosome and position
    #chrom, pos = line.split("\t")[0:2]
    chrom, pos = line.split()[0:2] #split on whitespace
    pos = int(filter(str.isalnum, pos)) #convert to integer
    #make sure the proper fastq is loaded
    if (chrom != curChrom):
        #load the entire fastq chrom
        faFile = os.path.join(faFol, chrom+".fa")        
        ref = loadFasta(faFile)
        #current
        curChrom = chrom
    #if the current position is a cpg, then print it
    #print ref[pos:(pos+2)]
    b2 = ref[pos:(pos+2)].lower()
    if b2 == "cg":
        #out.write("%s\t%s\n"%(line, b2))
        sys.stdout.write("%s\t%s\n"%(line, b2))
    
#out.close()



