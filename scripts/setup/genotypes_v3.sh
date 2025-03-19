
## Get all rat IDs for Adipose and Liver, which are in round 10.4 but not 10.5
# cat v2/{Adipose,Liver}/rat_ids* | sort | uniq | sed 's/_/-/' > geno/rfids_for_geno_10_4.txt
## Get all rat IDs from pre-QC and post-QC lists in v2 except Adipose and Liver, plus new tissues
# cat v2/{BLA,Brain,Eye,IL,LHb,NAcc,NAcc2,OFC,PL,PL2}/rat_ids* v3/{RMTg,NAcc3,PL3,pVTA}/rat_ids* | sort | uniq > geno/rfids_for_geno_10_5.txt
## Then on TSCC:
# awk 'FNR==NR { a[$0]; next } $0 in a' ../round10_4_rfids.txt rfids_for_gene_10_4.txt > rfids_in_round10_4.txt
# awk 'FNR==NR { a[$0]; next } $0 in a' ../round10_5_rfids.txt rfids_for_gene_10_5.txt > rfids_in_round10_5.txt

