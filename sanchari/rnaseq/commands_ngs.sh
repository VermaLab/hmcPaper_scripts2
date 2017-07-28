
#tss, tes, genebody, exon, cgi, enhancer, dhs
jdbam=JRN.AD0C99ABXX.lane_8_P0_I2.hg19.gsnap-rna-alignment.bam
hpnebam=HPNE.AC1FL7ACXX.lane_1_P0_I5.hg19.gsnap-rna-alignment.bam

#ngs.plot.r -P 4 -G hg19 -R genebody -C $jdbam -O jdplot_genebody -T jdRna -L 5000 
#ngs.plot.r -P 4 -G hg19 -R genebody -C $hpnebam -O hpneplot_genebody -T hpneRna -L 5000 

#
ngs.plot.r -P 4 -G hg19 -R tss -C $jdbam -O jdplot_rna_tss -T jdRna -L 5000 
ngs.plot.r -P 4 -G hg19 -R tes -C $jdbam -O jdplot_rna_tes -T jdRna -L 5000 
ngs.plot.r -P 4 -G hg19 -R exon -C $jdbam -O jdplot_rna_exon -T jdRna -L 5000 
ngs.plot.r -P 4 -G hg19 -R cgi -C $jdbam -O jdplot_rna_cgi -T jdRna -L 5000 
#
ngs.plot.r -P 4 -G hg19 -R tss -C $hpnebam -O hpneplot_rna_tss -T hpneRna -L 5000 
ngs.plot.r -P 4 -G hg19 -R tes -C $hpnebam -O hpneplot_rna_tes -T hpneRna -L 5000 
ngs.plot.r -P 4 -G hg19 -R exon -C $hpnebam -O hpneplot_rna_exon -T hpneRna -L 5000 
ngs.plot.r -P 4 -G hg19 -R cgi -C $hpnebam -O hpneplot_rna_cgi -T hpneRna -L 5000 

##run these two again
#ngs.plot.r -P 4 -G hg19 -R enhancer -C $jdbam -O jdplot_enhancer -T jdAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R dhs -C $jdbam -O jdplot_dhs -T jdAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R enhancer -C $hpnebam -O hpneplot_enhancer -T hpneAtac -L 5000 
#ngs.plot.r -P 4 -G hg19 -R dhs -C $hpnebam -O hpneplot_dhs -T hpneAtac -L 5000 

