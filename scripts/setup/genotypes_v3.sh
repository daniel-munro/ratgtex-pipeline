set -euxo pipefail

# This documents how the original dataset genotypes were processed uniformly.
# All genotypes are extracted from the round 10.5 VCF except the Adipose and
# Liver cohort, which are not in 10.5, so those are extracted from 10.4.

## Preprocessing:
## Get all rat IDs for Adipose and Liver, which are in round 10.4 but not 10.5
# cat v2/{Adipose,Liver}/rat_ids* | sort | uniq | sed 's/_/-/' > geno/original/rfids_for_geno_10_4.txt
## Get all rat IDs from pre-QC and post-QC lists in v2 except Adipose and Liver, plus new tissues
# cat v2/{BLA,Brain,Eye,IL,LHb,NAcc,NAcc2,OFC,PL,PL2}/rat_ids* v3/{RMTg,NAcc3,PL3,pVTA}/rat_ids.txt | sort | uniq > geno/original/rfids_for_geno_10_5.txt
## Then on TSCC:
# awk 'FNR==NR { a[$0]; next } $0 in a' ../round10_4_rfids.txt rfids_for_geno_10_4.txt > rfids_in_round10_4.txt
# awk 'FNR==NR { a[$0]; next } $0 in a' ../round10_5_rfids.txt rfids_for_geno_10_5.txt > rfids_in_round10_5.txt
# bcftools view /tscc/projects/ps-palmer/hs_rats/round10_4/genotypes/round10_4.vcf.gz --samples-file rfids_in_round10_4.txt | bcftools annotate -x ^INFO/AC,INFO/AN,^FORMAT/GT -Oz -o ratgtex_v3_round10_4.vcf.gz
# bcftools view /tscc/projects/ps-palmer/hs_rats/round10_5/genotypes/round10_5.vcf.gz --samples-file rfids_in_round10_5.txt | bcftools annotate -x ^INFO/AC,INFO/AN,^FORMAT/GT -Oz -o ratgtex_v3_round10_5.vcf.gz
## Then copy those to geno/original/

# Two steps are used to make sure REF alleles are all correct.
# 1. plink2 with '--ref-from-fa' swaps alleles and genotypes if ALT matches the true REF.
#    But if neither REF nor ALT match the true REF, it keeps the SNP as-is.
# 2. bcftools norm with '--check-ref x' removes those remaining mismatched SNPs.
# bcftools norm with '--check-ref s' will not work because if neither REF nor ALT match the true REF,
#    it keeps the SNP and just replaces the REF allele. But we shouldn't trust those genotypes.

mkdir -p geno/intermediate

### Round 10.4 genotypes ###
echo '*** Preparing round 10.4 genotypes...'
plink2 --vcf geno/original/ratgtex_v3_round10_4.vcf.gz \
    --chr 1-20 \
    --set-all-var-ids '@:#' \
    --fa ref/GCF_015227675.2_mRatBN7.2_genomic.chr.fa \
    --ref-from-fa force \
    --recode vcf \
    --output-chr chrM \
    --out geno/intermediate/ratgtex_v3_round10_4
bgzip geno/intermediate/ratgtex_v3_round10_4.vcf
bcftools norm geno/intermediate/ratgtex_v3_round10_4.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/GCF_015227675.2_mRatBN7.2_genomic.chr.fa \
    -Ou | bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    -O z -o geno/ratgtex_v3_round10_4.vcf.gz
tabix -f geno/ratgtex_v3_round10_4.vcf.gz

### Round 10.5 genotypes ###
echo '*** Preparing round 10.5 genotypes...'
plink2 --vcf geno/original/ratgtex_v3_round10_5.vcf.gz \
    --chr 1-20 \
    --set-all-var-ids '@:#' \
    --fa ref/GCF_015227675.2_mRatBN7.2_genomic.chr.fa \
    --ref-from-fa force \
    --recode vcf \
    --output-chr chrM \
    --out geno/intermediate/ratgtex_v3_round10_5
bgzip geno/intermediate/ratgtex_v3_round10_5.vcf
bcftools norm geno/intermediate/ratgtex_v3_round10_5.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/GCF_015227675.2_mRatBN7.2_genomic.chr.fa \
    -Ou | bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    -O z -o geno/ratgtex_v3_round10_5.vcf.gz
tabix -f geno/ratgtex_v3_round10_5.vcf.gz

### Get exon SNPs for all genotyped rats for sample mixup QC ###
zcat ref/exon_regions.tsv.gz | awk '$1 ~ /^chr([1-9]|1[0-9]|20)$/' > geno/intermediate/all_rats_exon_regions.tsv
## Then copy geno/intermediate/all_rats_exon_regions.tsv to TSCC and run:
# plink2 --vcf /tscc/projects/ps-palmer/hs_rats/round10_5/genotypes/round10_5.vcf.gz \
#     --chr 1-20 \
#     --extract bed1 all_rats_exon_regions.tsv \
#     --set-all-var-ids '@:#' \
#     --recode vcf \
#     --output-chr chrM \
#     --out tmp
# bgzip tmp.vcf
# bcftools annotate tmp.vcf.gz \
#     -x ^INFO/AC,INFO/AN,^FORMAT/GT \
#     -O z -o all_rats_exons.round10_5.vcf.gz
# rm tmp.vcf.gz
## Then copy all_rats_exons.round10_5.vcf.gz back to geno/all_rats_exons.vcf.gz and run:
tabix -f geno/all_rats_exons.vcf.gz
zcat geno/all_rats_exons.vcf.gz | grep -v '^#' | cut -f3 > geno/all_rats_exons.snps.txt

### Get all alleles ###
## -m none : no new multiallelics, output multiple records instead
## --force-samples : if the merged files contain duplicate samples names, proceed anyway
## awk '!_[$1]++' keeps only the first instance per SNP ID (in rare cases there are multiple)
bcftools merge \
    -m none \
    --force-samples \
    geno/ratgtex_v3_round10_4.vcf.gz \
    geno/ratgtex_v3_round10_5.vcf.gz \
    -Ou | bcftools view \
    --drop-genotypes \
    -Ov | grep -v '^#' | cut -f3-5 | awk '!_[$1]++' | gzip -c > geno/alleles.txt.gz

### Combine genotyping logs ###
cat geno/original/hs_rats_round10_4_n20031_20240624_genotype_log.csv <(cut -d',' -f1-36 geno/original/hs_rats_round10_5_n21864_20250227_genotype_log.csv | tail -n+2) > geno/genotyping_log.csv
