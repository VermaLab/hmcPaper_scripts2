options(stringsAsFactors=F)

COMMON_THRESH = 4
#combine the hpne and jd samples 

jd.files = paste0("hmcSites_JD_",c(1:22, "X", "Y"), "_cpg.txt")
hpne.files = paste0("hmcSites_HPNE_",c(1:22, "X", "Y"), "_cpg.txt")
out.files = paste0("hmcSites_BOTH_",c(1:22, "X", "Y"), "_cpg_common", COMMON_THRESH, ".txt")


for (i in 1:length(jd.files)){
    x.jd = read.table(jd.files[i])
    x.hpne = read.table(hpne.files[i])
    #head(x.jd)
    x.both = merge(x.jd, x.hpne, by=2)
    #head(x.both)

    #column names
    # chr
    # pos
    # #meth BS JD
    # #tot BS JD
    # #meth OXBS JD
    # #tot OXBS JD
    # odd ratio JD
    # fisher p-val JD
    # #meth BS HPNE
    # #tot BS HPNE
    # #meth OXBS HPNE
    # #tot OXBS HPNE
    # odd ratio HPNE
    # fisher p-val HPNE
    y = x.both[,c(2, 1, 3:8, 11:16)]
    head(y)
        #  V1.x    V2 V3.x V4.x V5.x V6.x V7.x V8.x V3.y V4.y V5.y V6.y V7.y V8.y
        #1 chr1 12264    3    3    0    1  Inf 0.25    1    1    3    3  NaN  1.0
        #2 chr1 12270    3    3    1    1  NaN 1.00    1    1    3    4  Inf  0.8
        #3 chr1 12277    3    3    1    1  NaN 1.00    0    1    4    4    0  1.0

    #the read counts for the four files
    col.ix = c(4, 6, 10, 12)
    #all four files must meet criteria
    row.ix = apply(y[,col.ix] >= COMMON_THRESH, MARGIN=1, all)

    #head(y[row.ix,])
    #dim(y[row.ix,])
    write.table(y[row.ix,], file=out.files[i], row.names=F, col.names=F, quote=F, sep="\t")
}














