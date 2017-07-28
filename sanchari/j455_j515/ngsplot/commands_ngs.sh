
#tss, tes, genebody, exon, cgi, enhancer, dhs
jdbam=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/jd_sorted_mdup_merged.bam
hpnebam=/home/kpradhan/Desktop/hpc_home/projects/sanchari/j455_j515/hpne_sorted_mdup.bam

#
#ngs.plot.r -P 4 -G hg19 -R tss -C $jdbam -O jdplot_tss -T jdAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R tes -C $jdbam -O jdplot_tes -T jdAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R genebody -C $jdbam -O jdplot_genebody -T jdAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R exon -C $jdbam -O jdplot_exon -T jdAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R cgi -C $jdbam -O jdplot_cgi -T jdAtac -L 5000 
#
#ngs.plot.r -P 4 -G hg19 -R tss -C $hpnebam -O hpneplot_tss -T hpneAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R tes -C $hpnebam -O hpneplot_tes -T hpneAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R genebody -C $hpnebam -O hpneplot_genebody -T hpneAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R exon -C $hpnebam -O hpneplot_exon -T hpneAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R cgi -C $hpnebam -O hpneplot_cgi -T hpneAtac -L 5000 

#run these two again
ngs.plot.r -P 4 -G hg19 -R enhancer -C $jdbam -O jdplot_enhancer -T jdAtac -L 5000 
ngs.plot.r -P 4 -G hg19 -R dhs -C $jdbam -O jdplot_dhs -T jdAtac -L 5000 
ngs.plot.r -P 4 -G hg19 -R enhancer -C $hpnebam -O hpneplot_enhancer -T hpneAtac -L 5000 
ngs.plot.r -P 4 -G hg19 -R dhs -C $hpnebam -O hpneplot_dhs -T hpneAtac -L 5000 

