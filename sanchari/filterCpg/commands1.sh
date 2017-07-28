##sort the files first
#sort -k1,1 -k2,2n /home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/commonSites/common4.txt > sorted_common4.txt
#                                                                                                                         
#sort  -k1,1 -k2,2n /home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/hpc_mcSites/common4_mcSites_HPNE.txt > sorted_common4_mcSites_HPNE.txt
#sort  -k1,1 -k2,2n /home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/hpc_mcSites/common4_mcSites_JD.txt > sorted_common4_mcSites_JD.txt
#                                                                                                                         
#sort  -k1,1 -k2,2n /home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/hmcSites/common4_hmcSites_HPNE.txt > sorted_common4_hmcSites_HPNE.txt
#sort  -k1,1 -k2,2n /home/kpradhan/Desktop/hpc_home/projects/sanchari/j432_j446/hmcSites/common4_hmcSites_JD.txt > sorted_common4_hmcSites_JD.txt




#
#time cat sorted_common4.txt | ./filterCpg.py >| common4cpg.txt
#
#time cat sorted_common4_mcSites_HPNE.txt | ./filterCpg.py >| common4cpg_mcSites_HPNE.txt
#time cat sorted_common4_mcSites_JD.txt | ./filterCpg.py >| common4cpg_mcSites_JD.txt
#
#time cat sorted_common4_hmcSites_HPNE.txt | ./filterCpg.py >| common4cpg_hmcSites_HPNE.txt
#time cat sorted_common4_hmcSites_JD.txt | ./filterCpg.py >| common4cpg_hmcSites_JD.txt


#get the diff sets
#sites in hpne but no jd

##mc
#comm -23 --check-order <(cut -f1-2 common4cpg_mcSites_JD.txt | sort) <(cut -f1-2 common4cpg_mcSites_HPNE.txt | sort) >| common4cpg_mcSites_JD-HPNE.txt
#comm -13 --check-order <(cut -f1-2 common4cpg_mcSites_JD.txt | sort) <(cut -f1-2 common4cpg_mcSites_HPNE.txt | sort) >| common4cpg_mcSites_HPNE-JD.txt
#
##hmc
#comm -23 --check-order <(cut -f1-2 common4cpg_hmcSites_JD.txt | sort) <(cut -f1-2 common4cpg_hmcSites_HPNE.txt | sort) >| common4cpg_hmcSites_JD-HPNE.txt
#comm -13 --check-order <(cut -f1-2 common4cpg_hmcSites_JD.txt | sort) <(cut -f1-2 common4cpg_hmcSites_HPNE.txt | sort) >| common4cpg_hmcSites_HPNE-JD.txt



#get #c's  #cpgs in common for JD
#get #meth's for c's and cpgs

#merge with common, and cpg
mc_sites=/home/kpradhan/mnt/hpc_home/projects/sanchari/j432_j446/hpc_mcSites/mcSites_JD_all.txt
common_cpg_sites=sorted_JDcommon4cpg.txt
common_c_sites=/home/kpradhan/mnt/hpc_home/projects/sanchari/j432_j446/commonSites/sorted_common4_JD.txt

ls $common_c_sites
ls $common_cpg_sites

head $mc_sites
head $common_c_sites
head $common_cpg_sites

cat $mc_sites | awk '{printf("%s\t%d#%s\n", $1,$2,$0)}' | head


cat $mc_sites | awk '{printf("%s\t%d#%s\n", $1,$2,$0)}' | head

cat $common_cpg_sites | cut -f1-2

cat $mc_sites | awk '{printf("%s\t%d#%s\n", $1,$2,$0)}' | head

#default sort doesn't work properly!
sort -t"#" -k1b,1 mc1.txt > mc2.txt
cat $common_cpg_sites | cut -f1-2 | sort >| cpg1.txt
cat $common_c_sites | cut -f1-2 | sort >| c1.txt

#cpg sites
join mc2.txt cpg1.txt -t'#' >| JdMcSites_JDcommon4Cpg.txt
join mc2.txt c1.txt -t'#' >| JdMcSites_JDcommon4C.txt

wc -l JdMcSites_JDcommon4Cpg.txt
wc -l JdMcSites_JDcommon4C.txt

join mc1.txt <(cat $common_cpg_sites | cut -f1-2 | sort) -t "#" >| JdMcSites_JDcommon4Cpg.txt

#c sites
join <(cat $mc_sites | awk '{printf("%s\t%d#%s\n", $1,$2,$0)}' | sort -t'#' -k1) <(cat $common_c_sites | cut -f1-2 | sort) -t "#"  >| JdMcSites_JDcommon4C.txt

