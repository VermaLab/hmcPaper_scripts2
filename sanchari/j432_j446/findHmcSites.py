#!/usr/bin/python

import gzip
import sys
from scipy.stats import fisher_exact

#ignore sites that have fewer than this many reads in both BS and OxBS
COMMON_THRESH = 4

#file are actually the unconverted counts, NOT converted counts
#   they show for each position
#   chr pos #meth #total
#  here, #meth are the number of unconverted reads!!

file1 = "conversionCounts/JD_BS_17.txt.gz"
file2 = "conversionCounts/JD_OxBs_17.txt.gz"
out = "hmcSites_py/hmcSites_HPNE_22.txt"
out = "hmcSites_py/test_22.txt"
out = "test1/test_17.txt"

file1 = sys.argv[1]
file2 = sys.argv[2]
out = sys.argv[3]

print (file1)
print (file2)
print (out)


fp1 = iter(gzip.open(file1, 'r'))
fp2 = iter(gzip.open(file2, "r"))
fp3 = open(out, "w")


x1 = next(fp1, "").rstrip("\n").split("\t")
x2 = next(fp2, "").rstrip("\n").split("\t")


while(True and len(x1) != 1 and len(x2) != 1):
    if (x1[0] != x2[0]):
        break
    #matching site, process and print
    if (int(x1[1]) == int(x2[1])):
        a = int(x1[2])
        b = int(x1[3]) - a
        c = int(x2[2])
        d = int(x2[3]) - c
        if (a+b >= COMMON_THRESH and c+b >= COMMON_THRESH):
            #r uses conditional MLE method to calc odds ratio
            #so it won't be exactly the same
            #odds, pval = fisher_exact([[a, b], [c, d]])
            #odds, pval = fisher_exact([[a, c], [b, d]], alternative="less")  #this one is incorrect!!
            #we're looking for high methylation in BS
            #and low methylation in OxBS
            odds, pval = fisher_exact([[a, c], [b, d]], alternative="greater")
            #process fisher test pvalue
            #
            #print result
            fp3.write("\t".join(map(str, (x1+x2[2:4] + [odds, pval])))+"\n")
            #if (odds < 1 and pval < 0.05):
            #    print odds, pval
            #    print x1
            #    print x2
            #    sigCount = sigCount + 1
            #
        #move both pointers up
        x1 = next(fp1, "").rstrip("\n").split("\t")
        x2 = next(fp2, "").rstrip("\n").split("\t")
    elif(int(x1[1]) < int(x2[1])):
        #move the bs file pointer up a spot
        x1 = next(fp1, "").rstrip("\n").split("\t")
    elif (int(x1[1]) > int(x2[1])):
        #move the oxbs file pointer up a spot
        x2 = next(fp2, "").rstrip("\n").split("\t")

fp1.close()
fp2.close()
fp3.close()
