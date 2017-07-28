#!/usr/bin/env Rscript

options(stringsAsFactors=F)

#command line arguments
args = commandArgs(trailingOnly=TRUE)
f.bs = args[1]
f.oxbs = args[2]
f.hmc = args[3]
#f.bs = "conversionCounts/HPNE_BS_22.txt.gz" 
#f.oxbs = "conversionCounts/HPNE_OxBS_22.txt.gz"
#f.hmc = "hmcSites_HPNE_22.txt"
#f.hmc = "test.txt"

print(f.bs)
print(f.oxbs)
print(f.hmc)

#BS
#c -> t
#5mc -> c
#5hmc -> c
#
#oxBS
#c -> t
#5mc -> c
#5hmc -> t

#we want to find the sites with hmc
#low conversion in BS
#high conversion in oxBS


#candidate sites should have 
#conversion rate in Bs  below 50%
#conversion rate in oxBs above 50%



#instead of loading everything in at once
#read and process line by line
con1 = pipe(paste0("zcat ", f.bs))
open(con1, "r")
con2 = pipe(paste0("zcat ", f.oxbs))
open(con2, "r")

#start the file pointer at the begining
x1 = read.table(textConnection(readLines(con1, n=1)))
x2 = read.table(textConnection(readLines(con2, n=1)))

#keep processing line by line until there's an error
try({while(T){
	#we should only be working with a single chrom
	if (x1[1,1] != x2[1,1]){
		break
	}

	#matching site, process and print
	if (x1[1,2] == x2[1,2]){
		x = merge(x1, x2, by.x=1, by.y=1)
		v =as.numeric(x[1, c("V3.x", "V4.x", "V3.y", "V4.y")])

		#if both files have at least 4 reads
		#and both have a 0 conversion rate
		#if (v[1] == 0 && v[2] >= 4 && v[3] == 0 && v[4] >= 4){ #this is incorrect
                #v1 is the # of unconverted reads of file1, v2 total file1
                #v3 is the # of unconverted reads of file2, v4 total file2
		if (v[1] == v[2] && v[2] >= 4 && v[3] == v[4] && v[4] >= 4){
			cat(paste(c(as.character(x[1,1:2]), v)), file=f.hmc, append=T)
			cat("\n", file=f.hmc, append=T)
		}
		#move both pointers up
		x1 = read.table(textConnection(readLines(con1, 1)))
		x2 = read.table(textConnection(readLines(con2, 1)))
	}else if(x1[1,2] < x2[1,2]){
		#move the bs file pointer up a spot
		x1 = read.table(textConnection(readLines(con1, 1)))
	}else if (x1[1,2] > x2[1,2]){
		#move the oxbs file pointer up a spot
		x2 = read.table(textConnection(readLines(con2, 1)))
	}
}}, silent=T)

close(con1)
close(con2)



