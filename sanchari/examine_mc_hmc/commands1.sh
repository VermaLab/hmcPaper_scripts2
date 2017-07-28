#prepare the oxbs files for a specific region.

#cut out the specific regions belonging to the gene.

#plus promoters
#chr19:15235519-15332545

#or just read each chromosome in, and use R to get the region chunk

#lets try lookin at just the cpg sites
cat hpne_bs.bdg | python ../filterCpg/filterCpg.py > hpne_bs_cpg.bdg
cat hpne_oxbs.bdg | python ../filterCpg/filterCpg.py > hpne_oxbs_cpg.bdg
cat jd_bs.bdg | python ../filterCpg/filterCpg.py > jd_bs_cpg.bdg
cat jd_oxbs.bdg | python ../filterCpg/filterCpg.py > jd_oxbs_cpg.bdg

#filter all but the cpg sites and save into new files
ls *.bdg | grep -v -e "-cpg." | sed 's/\(.*\)\.bdg/cat & \| python \.\.\/filterCpg\/filterCpg\.py >| \1_cpg.bdg/'

#merge all the cpg files by sample type
cat *_hpne_oxbs_cpg* >| convPerc_hpne_oxbs_cpg.bdg
cat *_hpne_bs_cpg* >| convPerc_hpne_bs_cpg.bdg
cat *_jd_oxbs_cpg* >| convPerc_jd_oxbs_cpg.bdg
cat *_jd_bs_cpg* >| convPerc_jd_bs_cpg.bdg








