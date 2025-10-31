set -e

###################
## Set up for QC ##
###################

while read tissue; do
    mkdir -p v4/${tissue}
    ## Some of these won't work and must be done manually
    cp v3/${tissue}/rat_ids.txt v4/${tissue}/
    cp v3/${tissue}/fastq_map.txt v4/${tissue}/
done < tissues.dup.txt

mkdir -p geno/original
cat v4/*/rat_ids.txt | sort | uniq > geno/original/rfids_for_geno_11_1.txt

## Transfer round11_1_vcf_ids.txt from TSCC. These IDs have extra info around RFIDs and must be parsed.

## Make lists of RFIDs and corresponding VCF IDs to extract genotypes
python scripts/setup/genotype_ids.py \
    --rfids geno/original/rfids_for_geno_11_1.txt \
    --vcf-ids geno/original/round11_1_vcf_ids.txt \
    --out-tsv geno/original/id_map_to_extract_11_1.tsv

## Transfer output to TSCC and run the commands documented in comments in genotypes_v4.sh

#######################
## Set up for Pantry ##
#######################

while read tissue; do
    rsync -av ~/tools/Pantry/phenotyping/ v4/$tissue/phenos --exclude input --exclude intermediate --exclude output --exclude .snakemake
    cp scripts/config_phenos.yml v4/$tissue/phenos/config.yml
    # (Then update read length in config.yml)
    mkdir -p v4/$tissue/phenos/intermediate/bam
done < tissues.dup.txt

## Link FASTQ and BAM files

ln -s fastq/BLA_NAcc2_PL2 v4/BLA/fastq
ln -s fastq/Brain v4/Brain/fastq
ln -s fastq/Eye v4/Eye/fastq
ln -s fastq/IC_NAcc4_pVTA2 v4/IC/fastq
ln -s fastq/IL_LHb_NAcc_OFC_PL v4/IL/fastq
ln -s fastq/IL_LHb_NAcc_OFC_PL v4/LHb/fastq
ln -s fastq/IL_LHb_NAcc_OFC_PL v4/NAcc1/fastq
ln -s fastq/BLA_NAcc2_PL2 v4/NAcc2/fastq
ln -s fastq/NAcc3_PL3_pVTA v4/NAcc3/fastq
ln -s fastq/IC_NAcc4_pVTA2 v4/NAcc4/fastq
ln -s fastq/IL_LHb_NAcc_OFC_PL v4/OFC/fastq
ln -s fastq/IL_LHb_NAcc_OFC_PL v4/PL1/fastq
ln -s fastq/BLA_NAcc2_PL2 v4/PL2/fastq
ln -s fastq/NAcc3_PL3_pVTA v4/PL3/fastq
ln -s fastq/NAcc3_PL3_pVTA1 v4/pVTA1/fastq
ln -s fastq/IC_NAcc4_pVTA2 v4/pVTA2/fastq
ln -s fastq/RMTg v4/RMTg/fastq

while read tissue; do
    while read sample; do
        ln -s ../../../star_out/$sample.Aligned.sortedByCoord.out.bam v4/$tissue/phenos/intermediate/bam/$sample.bam
        ln -s ../../../star_out/$sample.Aligned.sortedByCoord.out.bam.bai v4/$tissue/phenos/intermediate/bam/$sample.bam.bai
    done < v4/$tissue/rat_ids.txt
done < tissues.dup.txt

