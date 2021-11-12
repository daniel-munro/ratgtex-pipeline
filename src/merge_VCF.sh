set -euxo pipefail

# This documents how the original dataset genotypes were merged into geno/ratgtex.vcf.gz.
# Two steps are used to make sure REF alleles are all correct.
# 1. plink2 with '--ref-from-fa' swaps alleles and genotypes if ALT matches the true REF.
#    But if neither REF nor ALT match the true REF, it keeps the SNP as-is.
# 2. bcftools norm with '--check-ref x' removes those remaining mismatched SNPs.
# bcftools norm with '--check-ref s' will not work because if neither REF nor ALT match the true REF,
#    it keeps the SNP and just replaces the REF allele. But we shouldn't trust those genotypes.

BRAIN_DIR=~/br/data/genotype
EYE_DIR=~/eye/data/genotype
ADIPOSE_DIR=~/bulk/fl/Imputed_Geno

mkdir -p geno/intermediate

### Brain ###
echo '*** Preparing Brain genotypes...'
plink2 --vcf $BRAIN_DIR/P50.rnaseq.88.unpruned.vcf.gz \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno/intermediate/brain1
bgzip geno/intermediate/brain1.vcf
bcftools norm geno/intermediate/brain1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -O z -o geno/intermediate/brain.vcf.gz
tabix -f geno/intermediate/brain.vcf.gz

### Eye ###
echo '*** Preparing Eye genotypes...'
plink2 --vcf $EYE_DIR/eyes.vcf.gz \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno/intermediate/eyes1
bgzip geno/intermediate/eyes1.vcf
bcftools norm geno/intermediate/eyes1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -O z -o geno/intermediate/eyes.vcf.gz
tabix -f geno/intermediate/eyes.vcf.gz

### Adipose ###
echo '*** Preparing Adipose genotypes...'
plink2 --bfile $ADIPOSE_DIR/RNAed.chrall \
    --recode vcf id-paste=iid \
    --set-all-var-ids 'chr@:#' \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --out geno/intermediate/adipose1
bgzip geno/intermediate/adipose1.vcf
bcftools norm geno/intermediate/adipose1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -O z -o geno/intermediate/adipose.vcf.gz
tabix -f geno/intermediate/adipose.vcf.gz

### Merge ###
echo '*** Merging across datasets...'
bcftools merge \
    geno/intermediate/brain.vcf.gz \
    geno/intermediate/eyes.vcf.gz \
    geno/intermediate/adipose.vcf.gz \
    -O z -o geno/ratgtex_tmp.vcf.gz
# Reheader for correct contig lengths and to remove confusing extra info:
cp geno/header.txt geno/header_tmp.txt
bcftools view -h geno/ratgtex_tmp.vcf.gz | grep '^#CHROM' >> geno/header_tmp.txt
bcftools reheader -h geno/header_tmp.txt geno/ratgtex_tmp.vcf.gz -o geno/ratgtex.vcf.gz
rm geno/ratgtex_tmp.vcf.gz geno/header_tmp.txt
tabix -f geno/ratgtex.vcf.gz

### Troubleshooting

# #Get SNPs with multiple REF alleles (i.e. wrong allele in at least one dataset)
# cat \
#     <(zcat geno/intermediate/brain.vcf.gz | grep -v '^#' | cut -f3,4) \
#     <(zcat geno/intermediate/eyes.vcf.gz | grep -v '^#' | cut -f3,4) \
#     <(zcat geno/intermediate/adipose.vcf.gz | grep -v '^#' | cut -f3,4) \
#     | sort | uniq | cut -f1 | sort | uniq -d > dup.txt

# # Get SNP 'regions' to get true REF alleles using samtools faidx:
# zcat geno/intermediate/adipose.vcf.gz | grep -v '^#' | cut -f1,2 \
#     | awk '{print $1 ":" $2 "-" $2}' | head > geno/intermediate/regions.adipose.txt

# # Get true REF allele for each of those SNPs
# samtools faidx ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa -r geno/intermediate/regions.adipose.txt
