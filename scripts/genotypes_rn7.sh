set -euxo pipefail

# This documents how the original dataset genotypes were processed uniformly.
# Two steps are used to make sure REF alleles are all correct.
# 1. plink2 with '--ref-from-fa' swaps alleles and genotypes if ALT matches the true REF.
#    But if neither REF nor ALT match the true REF, it keeps the SNP as-is.
# 2. bcftools norm with '--check-ref x' removes those remaining mismatched SNPs.
# bcftools norm with '--check-ref s' will not work because if neither REF nor ALT match the true REF,
#    it keeps the SNP and just replaces the REF allele. But we shouldn't trust those genotypes.

BRAIN_REGION_DIR=~/bulk/br/geno_rn7
ROUND10_VCF=~/ratgtex/geno_rn7/original/Heterogenous-stock_n15552_02222023_stitch2_QC_Sex_Het_pass_n14505.vcf.gz

mkdir -p geno_rn7/intermediate

### IL_LHb_NAcc_OFC_PL ###
echo '*** Preparing IL_LHb_NAcc_OFC_PL genotypes...'
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
tabix -f geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.2.vcf.gz
## Rescue one sample that was contaminated in the 30x genotypes and removed
## plus two samples that were mislabeled but matched genotypes in round 10
## (Do this plink step before merging to fix the chrom names being different formats)
plink2 --vcf $ROUND10_VCF \
    --keep geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.round10.txt \
    --set-all-var-ids 'chr@:#' \
    --fa ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.round10
bgzip geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.round10.vcf
tabix -f geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.round10.vcf.gz
bcftools merge \
    geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.2.vcf.gz \
    geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.round10.vcf.gz \
    -O z -o geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.3.vcf.gz
bcftools norm geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.3.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    -Ou | bcftools annotate \
    -x ^INFO/AC,INFO/AN,^FORMAT/GT \
    -O z -o geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz
tabix -f geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz

### Eye ###
echo '*** Preparing Eye genotypes...'
plink2 --vcf $ROUND10_VCF \
    --set-all-var-ids 'chr@:#' \
    --fa ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    --ref-from-fa force \
    --keep geno_rn7/intermediate/Eye_ids.txt \
    --recode vcf \
    --out geno_rn7/intermediate/Eye.1
bgzip geno_rn7/intermediate/Eye.1.vcf
bcftools norm geno_rn7/intermediate/Eye.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    -Ou | bcftools annotate \
    -x ^INFO/AC,INFO/AN,^FORMAT/GT \
    -O z -o geno_rn7/intermediate/Eye.vcf.gz
tabix -f geno_rn7/intermediate/Eye.vcf.gz

### Brain ###
echo '*** Preparing whole brain genotypes...'
plink2 --vcf $ROUND10_VCF \
    --set-all-var-ids 'chr@:#' \
    --fa ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    --ref-from-fa force \
    --keep geno_rn7/intermediate/Brain_ids.txt \
    --recode vcf \
    --out geno_rn7/intermediate/Brain.1
bgzip geno_rn7/intermediate/Brain.1.vcf
bcftools norm geno_rn7/intermediate/Brain.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    -Ou | bcftools annotate \
    -x ^INFO/AC,INFO/AN,^FORMAT/GT \
    -O z -o geno_rn7/intermediate/Brain.vcf.gz
tabix -f geno_rn7/intermediate/Brain.vcf.gz

### BLA_NAcc2_PL2 ###
echo '*** Preparing BLA_NAcc2_PL2 genotypes...'
plink2 --vcf $ROUND10_VCF \
    --set-all-var-ids 'chr@:#' \
    --fa ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    --ref-from-fa force \
    --keep geno_rn7/intermediate/BLA_NAcc2_PL2_ids.txt \
    --recode vcf \
    --out geno_rn7/intermediate/BLA_NAcc2_PL2.1
bgzip geno_rn7/intermediate/BLA_NAcc2_PL2.1.vcf
bcftools norm geno_rn7/intermediate/BLA_NAcc2_PL2.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa \
    -Ou | bcftools annotate \
    -x ^INFO/AC,INFO/AN,^FORMAT/GT \
    -O z -o geno_rn7/intermediate/BLA_NAcc2_PL2.vcf.gz
tabix -f geno_rn7/intermediate/BLA_NAcc2_PL2.vcf.gz

### Final processing ###
# for DSET in IL_LHb_NAcc_OFC_PL Eye Adipose_Liver Brain BLA_NAcc2_PL2; do
for DSET in IL_LHb_NAcc_OFC_PL; do
    echo "*** Final processing: $DSET..."
    bcftools view geno_rn7/intermediate/$DSET.vcf.gz \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        --targets "^Y,MT" \
        -O z -o geno_rn7/$DSET.tmp.vcf.gz
    # Reheader to remove confusing extra info:
    cp geno_rn7/header.txt geno_rn7/header.tmp.txt
    bcftools view -h geno_rn7/$DSET.tmp.vcf.gz | grep '^#CHROM' >> geno_rn7/header.tmp.txt
    bcftools reheader -h geno_rn7/header.tmp.txt geno_rn7/$DSET.tmp.vcf.gz -o geno_rn7/$DSET.vcf.gz
    rm geno_rn7/$DSET.tmp.vcf.gz geno_rn7/header.tmp.txt
    tabix -f geno_rn7/$DSET.vcf.gz
done

### Get all alleles ###
## -m none
## --force-samples
## awk '!_[$1]++' keeps only the first instance per SNP ID (in rare cases there are multiple)
bcftools merge \
    -m none \
    --force-samples \
    geno_rn7/intermediate/BLA_NAcc2_PL2.vcf.gz \
    geno_rn7/intermediate/Brain.vcf.gz \
    geno_rn7/intermediate/Eye.vcf.gz \
    geno_rn7/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz \
    -Ou | bcftools view \
    --drop-genotypes \
    -Ov | grep -v '^#' | cut -f3-5 | awk '!_[$1]++' | gzip -c > geno_rn7/alleles.txt.gz

### Get exon SNPs for all genotyped rats for sample mixup QC
zcat ref_rn7/exon_regions.tsv.gz | head -239746 > geno_rn7/intermediate/all_rats_exon_regions.tsv
plink2 --vcf $ROUND10_VCF \
    --extract bed1 geno_rn7/intermediate/all_rats_exon_regions.tsv \
    --set-all-var-ids 'chr@:#' \
    --recode vcf \
    --out geno_rn7/intermediate/all_rats_exons
bgzip geno_rn7/intermediate/all_rats_exons.vcf
bcftools annotate geno_rn7/intermediate/all_rats_exons.vcf.gz \
    -x ^INFO/AC,INFO/AN,^FORMAT/GT \
    -O z -o geno_rn7/all_rats_exons.vcf.gz
tabix -f geno_rn7/all_rats_exons.vcf.gz
