set -euxo pipefail

# This documents how the original dataset genotypes were processed uniformly.
# Two steps are used to make sure REF alleles are all correct.
# 1. plink2 with '--ref-from-fa' swaps alleles and genotypes if ALT matches the true REF.
#    But if neither REF nor ALT match the true REF, it keeps the SNP as-is.
# 2. bcftools norm with '--check-ref x' removes those remaining mismatched SNPs.
# bcftools norm with '--check-ref s' will not work because if neither REF nor ALT match the true REF,
#    it keeps the SNP and just replaces the REF allele. But we shouldn't trust those genotypes.

BRAIN_REGION_DIR=~/br/data/genotype
# EYE_DIR=~/eye/data/genotype
ADIPOSE_LIVER_DIR=~/bulk/fl/Imputed_Geno
WHOLEBR_DIR=~/wb/data/genotype/TWAS_Whole_Brain
MITCHELL_DIR=~/sm/data/genotype/Mitchell_Hitzemann

mkdir -p geno/intermediate

### IL_LHb_NAcc_OFC_PL ###
echo '*** Preparing Brain region genotypes...'
plink2 --vcf $BRAIN_REGION_DIR/P50.rnaseq.88.unpruned.vcf.gz \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno/intermediate/IL_LHb_NAcc_OFC_PL.1
bgzip geno/intermediate/IL_LHb_NAcc_OFC_PL.1.vcf
bcftools norm geno/intermediate/IL_LHb_NAcc_OFC_PL.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -O z -o geno/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz
tabix -f geno/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz

### Eye ###
echo '*** Preparing Eye genotypes...'
### This part was run in another script:
## for chr in {1..20}; do
##     echo $chr
##     rsync -av tscc:/projects/ps-palmer/apurva/phased_genotypes/chr${chr}/chr${chr}.P50.6147.unpruned.vcf .
##     bgzip chr${chr}.P50.6147.unpruned.vcf
##     bcftools view \
##         chr${chr}.P50.6147.unpruned.vcf.gz \
##         --samples-file ~/ratgtex/Eye/rat_ids.txt \
##         -Oz -o tmp_eye.chr${chr}.vcf.gz
## done
## # Then combine:
## bcftools concat \
##     tmp_eye.chr{1..20}.vcf.gz \
##     -Ou | bcftools annotate \
##     --rename-chrs chrs.txt \
##     -Oz -o ../intermediate/Eye.1.vcf.gz
### Now for actual merge:
plink2 --vcf geno/intermediate/Eye.1.vcf.gz \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno/intermediate/Eye.2
bgzip geno/intermediate/Eye.2.vcf
bcftools norm geno/intermediate/Eye.2.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -O z -o geno/intermediate/Eye.vcf.gz
tabix -f geno/intermediate/Eye.vcf.gz

### Adipose_Liver ###
echo '*** Preparing Adipose genotypes...'
plink2 --bfile $ADIPOSE_DIR/RNAed.chrall \
    --recode vcf id-paste=iid \
    --set-all-var-ids 'chr@:#' \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --out geno/intermediate/Adipose_Liver.1
bgzip geno/intermediate/Adipose_Liver.1.vcf
bcftools norm geno/intermediate/Adipose_Liver.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -O z -o geno/intermediate/Adipose_Liver.vcf.gz
tabix -f geno/intermediate/Adipose_Liver.vcf.gz

### Brain ###
echo '*** Preparing whole brain genotypes...'
plink2 --vcf $WHOLEBR_DIR/Heterogenous-stock_n4140_11152021_stitch1_NO_QC_TWAS_Whole_Brain_n322.vcf.gz \
    --set-all-var-ids 'chr@:#' \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno/intermediate/Brain.1
bgzip geno/intermediate/Brain.1.vcf
bcftools norm geno/intermediate/Brain.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -Ou | bcftools annotate \
    -x INFO/EAF,INFO/INFO_SCORE,INFO/HWE,INFO/ERC,INFO/EAC,INFO/PAF,INFO/REF_PANEL \
    -O z -o geno/intermediate/Brain.vcf.gz
tabix -f geno/intermediate/Brain.vcf.gz

### BLA_NAcc2_PL2 ###
echo '*** Preparing u01_suzanne_mitchell genotypes...'
plink2 --vcf $MITCHELL_DIR/Heterogenous-stock_n4140_11152021_stitch1_QC_sex_missing_het_pass_Mitchell_Hitzemann_n191.vcf.gz \
    --set-all-var-ids 'chr@:#' \
    --fa ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    --ref-from-fa force \
    --recode vcf \
    --out geno/intermediate/BLA_NAcc2_PL2.1
bgzip geno/intermediate/BLA_NAcc2_PL2.1.vcf
bcftools norm geno/intermediate/BLA_NAcc2_PL2.1.vcf.gz \
    --rm-dup snps \
    --check-ref x \
    --fasta-ref ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa \
    -Ou | bcftools annotate \
    -x INFO/EAF,INFO/INFO_SCORE,INFO/HWE,INFO/ERC,INFO/EAC,INFO/PAF,INFO/REF_PANEL \
    -O z -o geno/intermediate/BLA_NAcc2_PL2.vcf.gz
tabix -f geno/intermediate/BLA_NAcc2_PL2.vcf.gz

### Final processing ###
for DSET in IL_LHb_NAcc_OFC_PL Eye Adipose_Liver Brain BLA_NAcc2_PL2; do
    echo "*** Final processing: $DSET..."
    bcftools view geno/intermediate/$DSET.vcf.gz \
        --min-alleles 2 \
        --max-alleles 2 \
        --types snps \
        -O z -o geno/$DSET.tmp.vcf.gz
    # Reheader for correct contig lengths and to remove confusing extra info:
    cp geno/header.txt geno/header.tmp.txt
    bcftools view -h geno/$DSET.tmp.vcf.gz | grep '^#CHROM' >> geno/header.tmp.txt
    bcftools reheader -h geno/header.tmp.txt geno/$DSET.tmp.vcf.gz -o geno/$DSET.vcf.gz
    rm geno/$DSET.tmp.vcf.gz geno/header.tmp.txt
    tabix -f geno/$DSET.vcf.gz
done

### Get all alleles ###
## -m none
## --force-samples
## awk '!_[$1]++' keeps only the first instance per SNP ID (in rare cases there are multiple)
bcftools merge \
    -m none \
    --force-samples \
    geno/intermediate/Adipose_Liver.vcf.gz \
    geno/intermediate/BLA_NAcc2_PL2.vcf.gz \
    geno/intermediate/Brain.vcf.gz \
    geno/intermediate/Eye.vcf.gz \
    geno/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz \
    -Ou | bcftools view \
    --drop-genotypes \
    -Ov | grep -v '^#' | cut -f3-5 | awk '!_[$1]++' | gzip -c > geno/alleles.txt.gz

### Troubleshooting

# #Get SNPs with multiple REF alleles (i.e. wrong allele in at least one dataset)
# cat \
#     <(zcat geno/intermediate/IL_LHb_NAcc_OFC_PL.vcf.gz | grep -v '^#' | cut -f3,4) \
#     <(zcat geno/intermediate/Eye.vcf.gz | grep -v '^#' | cut -f3,4) \
#     <(zcat geno/intermediate/adipose.vcf.gz | grep -v '^#' | cut -f3,4) \
#     | sort | uniq | cut -f1 | sort | uniq -d > dup.txt

# # Get SNP 'regions' to get true REF alleles using samtools faidx:
# zcat geno/intermediate/adipose.vcf.gz | grep -v '^#' | cut -f1,2 \
#     | awk '{print $1 ":" $2 "-" $2}' | head > geno/intermediate/regions.adipose.txt

# # Get true REF allele for each of those SNPs
# samtools faidx ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa -r geno/intermediate/regions.adipose.txt


# ### Coding SNPs for all available HS rats
# # These are used to check for a genotype match for RNA samples still missing genotypes after mixup detection.
# # Run for each chromosome, e.g. as a bunch of jobs:
# chr=$1
# bcftools annotate \
#     chr${chr}.P50.6147.unpruned.vcf.gz \
#     --rename-chrs chrs.txt \
#     -Oz -o tmp1.chr${chr}.vcf.gz
# tabix -f tmp1.chr${chr}.vcf.gz
# bcftools view tmp1.chr${chr}.vcf.gz \
#     --regions-file ~/ratgtex/ref/exon_regions.tsv.gz \
#     -Oz -o tmp.chr${chr}.vcf.gz
# rm tmp1.chr${chr}.vcf.gz*
# # Then combine:
# bcftools concat \
#     tmp.chr{1..20}.vcf.gz \
#     -Oz -o all_rats_exons.vcf.gz
