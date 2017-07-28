
#get the sites in JD but not in HPNE
comm -23 --check-order <(cut -f1-2 common4_hmcSites_JD.txt | sort) <(cut -f1-2 common4_hmcSites_HPNE.txt | sort) >| common4_hmcSites_JD-HPNE.txt
comm -13 --check-order <(cut -f1-2 common4_hmcSites_JD.txt | sort) <(cut -f1-2 common4_hmcSites_HPNE.txt | sort) >| common4_hmcSites_HPNE-JD.txt




