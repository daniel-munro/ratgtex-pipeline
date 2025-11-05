set -euxo pipefail

# This documents how the original dataset genotypes were processed uniformly.
# All genotypes are extracted from the round 11.1 VCF.

## On TSCC:
# bcftools view \
#     /tscc/projects/ps-palmer/hs_rats/round11_1/output/08_sample_filter/round11_1_allchrom_allsexes_snp_and_sample_filtered.vcf.gz \
#     --samples-file <(cut -f1 id_map_to_extract_11_1.tsv) \
# | bcftools reheader -s id_map_to_extract_11_1.tsv \
# | bcftools annotate -x ^INFO/AC,INFO/AN,^FORMAT/GT \
#     -Oz -o ratgtex_v4_round11_1.vcf.gz
## Then copy that to geno/original/

# Two steps are used to make sure REF alleles are all correct.
# 1. plink2 with '--ref-from-fa' swaps alleles and genotypes if ALT matches the true REF.
#    But if neither REF nor ALT match the true REF, it keeps the SNP as-is.
# 2. bcftools norm with '--check-ref x' removes those remaining mismatched SNPs.
# bcftools norm with '--check-ref s' will not work because if neither REF nor ALT match the true REF,
#    it keeps the SNP and just replaces the REF allele. But we shouldn't trust those genotypes.

mkdir -p geno/intermediate

### Round 11.1 genotypes ###
echo '*** Preparing round 11.1 genotypes...'
plink2 --vcf geno/original/ratgtex_v4_round11_1.vcf.gz \
    --chr 1-20 \
    --set-all-var-ids '@:#' \
    --fa ref/GCF_036323735.1_GRCr8_genomic.chr.fa \
    --ref-from-fa force \
    --recode vcf \
    --output-chr chrM \
    --out geno/intermediate/ratgtex_v4_round11_1
bgzip geno/intermediate/ratgtex_v4_round11_1.vcf
bcftools norm geno/intermediate/ratgtex_v4_round11_1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/GCF_036323735.1_GRCr8_genomic.chr.fa \
    -Ou | bcftools view \
    --min-alleles 2 \
    --max-alleles 2 \
    -O z -o geno/ratgtex_v4_round11_1.vcf.gz
tabix -f geno/ratgtex_v4_round11_1.vcf.gz

### Get exon SNPs for all genotyped rats for sample mixup QC ###
zcat ref/exon_regions.tsv.gz | awk '$1 ~ /^chr([1-9]|1[0-9]|20)$/' > geno/intermediate/all_rats_exon_regions.tsv
## Then copy geno/intermediate/all_rats_exon_regions.tsv to TSCC and run:
# plink2 --vcf /tscc/projects/ps-palmer/hs_rats/round11_1/output/08_sample_filter/round11_1_allchrom_allsexes_snp_and_sample_filtered.vcf.gz \
#     --chr 1-20 \
#     --extract bed1 all_rats_exon_regions.tsv \
#     --set-all-var-ids '@:#' \
#     --recode vcf \
#     --output-chr chrM \
#     --out tmp
# bgzip tmp.vcf
# bcftools annotate tmp.vcf.gz \
#     -x ^INFO/AC,INFO/AN,^FORMAT/GT \
#     -O z -o all_rats_exons.round11_1.vcf.gz
# rm tmp.vcf.gz
## Then copy all_rats_exons.round11_1.vcf.gz back to geno/all_rats_exons.vcf.gz and run:
tabix -f geno/all_rats_exons.vcf.gz
zcat geno/all_rats_exons.vcf.gz | grep -v '^#' | cut -f3 > geno/all_rats_exons.snps.txt

### Get all alleles ###
## awk '!_[$1]++' keeps only the first instance per SNP ID (in rare cases there are multiple)
bcftools view \
    geno/ratgtex_v4_round11_1.vcf.gz \
    --drop-genotypes \
    -Ov | grep -v '^#' | cut -f3-5 | awk '!_[$1]++' | gzip -c > geno/alleles.txt.gz
