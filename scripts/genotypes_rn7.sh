set -euxo pipefail

# This documents how the original dataset genotypes were processed uniformly.
# Two steps are used to make sure REF alleles are all correct.
# 1. plink2 with '--ref-from-fa' swaps alleles and genotypes if ALT matches the true REF.
#    But if neither REF nor ALT match the true REF, it keeps the SNP as-is.
# 2. bcftools norm with '--check-ref x' removes those remaining mismatched SNPs.
# bcftools norm with '--check-ref s' will not work because if neither REF nor ALT match the true REF,
#    it keeps the SNP and just replaces the REF allele. But we shouldn't trust those genotypes.

BRAIN_REGION_DIR=~/bulk/br/geno_rn7

mkdir -p geno_rn7/intermediate

### IL_LHb_NAcc_OFC_PL ###
bcftools annotate $BRAIN_REGION_DIR/88_outbred_HS_mRatBN7_2.vcf.gz \
    -x ^INFO/AC,INFO/AN,^FORMAT/GT \
    --rename-chrs geno_rn7/chroms.tsv \
    -O z -o geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.1.vcf.gz
plink2 --vcf geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.1.vcf.gz \
    --set-all-var-ids 'chr@:#' \
    --fa ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.2
bgzip geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.2.vcf
bcftools norm geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.2.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    -O z -o geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz
tabix -f geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz

### Final processing ###
# for DSET in IL_LHb_NAcc_OFC_PL Eye Adipose_Liver Brain BLA_NAcc2_PL2; do
for DSET in IL_LHb_NAcc_OFC_PL; do
    echo "*** Final processing: $DSET..."
    bcftools view geno_rn7/intermediate/$DSET.vcf.gz \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        -O z -o geno_rn7/$DSET.tmp.vcf.gz
    # Reheader to remove confusing extra info:
    cp geno_rn7/header.txt geno_rn7/header.tmp.txt
    bcftools view -h geno_rn7/$DSET.tmp.vcf.gz | grep '^#CHROM' >> geno_rn7/header.tmp.txt
    bcftools reheader -h geno_rn7/header.tmp.txt geno_rn7/$DSET.tmp.vcf.gz -o geno_rn7/$DSET.vcf.gz
    rm geno_rn7/$DSET.tmp.vcf.gz geno_rn7/header.tmp.txt
    tabix -f geno_rn7/$DSET.vcf.gz
done

# ### Get all alleles ###
# ## -m none
# ## --force-samples
# ## awk '!_[$1]++' keeps only the first instance per SNP ID (in rare cases there are multiple)
# bcftools merge \
#     -m none \
#     --force-samples \
#     geno_rn7/intermediate/Adipose_Liver.vcf.gz \
#     geno_rn7/intermediate/BLA_NAcc2_PL2.vcf.gz \
#     geno_rn7/intermediate/Brain.vcf.gz \
#     geno_rn7/intermediate/Eye.vcf.gz \
#     geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz \
#     -Ou | bcftools view \
#     --drop-genotypes \
#     -Ov | grep -v '^#' | cut -f3-5 | awk '!_[$1]++' | gzip -c > geno_rn7/alleles.txt.gz
