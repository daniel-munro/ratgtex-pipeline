set -e

while read tissue; do
    mkdir -p v4/${tissue}
    # Some of these won't work and must be done manually
    cp v3/${tissue}/rat_ids.txt v4/${tissue}/
    cp v3/${tissue}/fastq_map.txt v4/${tissue}/
done < tissues.dup.txt

mkdir -p geno/original
## Include 35 rats missing from last round of genotypes
cat v4/*/rat_ids.txt \
    <(cut -f1 ../p50-hao-chen-2020/data/p50_hao_chen_2020_rnaseq_ids.txt) \
    | sort | uniq > geno/original/rfids_for_geno_11_1.txt
