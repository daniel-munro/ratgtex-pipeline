set -e

mkdir -p v4/{NAcc,PL,pVTA}
cat v4/NAcc{1,2,3,4}/rat_ids.txt | sort | uniq > v4/NAcc/rat_ids.txt
cat v4/PL{1,2,3}/rat_ids.txt | sort | uniq > v4/PL/rat_ids.txt
cat v4/pVTA{1,2}/rat_ids.txt | sort | uniq > v4/pVTA/rat_ids.txt

## fastq_map.txt are created just to provide for download
## Prepend sub-tissue ID to each FASTQ path
awk '{for(i=1;i<NF;i++) printf "NAcc1/%s\t", $i; print $NF}' v4/NAcc1/fastq_map.txt > v4/NAcc/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "NAcc2/%s\t", $i; print $NF}' v4/NAcc2/fastq_map.txt >> v4/NAcc/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "NAcc3/%s\t", $i; print $NF}' v4/NAcc3/fastq_map.txt >> v4/NAcc/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "NAcc4/%s\t", $i; print $NF}' v4/NAcc4/fastq_map.txt >> v4/NAcc/fastq_map.txt

awk '{for(i=1;i<NF;i++) printf "PL1/%s\t", $i; print $NF}' v4/PL1/fastq_map.txt > v4/PL/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "PL2/%s\t", $i; print $NF}' v4/PL2/fastq_map.txt >> v4/PL/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "PL3/%s\t", $i; print $NF}' v4/PL3/fastq_map.txt >> v4/PL/fastq_map.txt

awk '{for(i=1;i<NF;i++) printf "pVTA1/%s\t", $i; print $NF}' v4/pVTA1/fastq_map.txt > v4/pVTA/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "pVTA2/%s\t", $i; print $NF}' v4/pVTA2/fastq_map.txt >> v4/pVTA/fastq_map.txt

for tissue in NAcc PL pVTA; do
    rsync -av ~/tools/Pantry/phenotyping/ v4/$tissue/phenos --exclude input --exclude intermediate --exclude output --exclude .snakemake
    cp scripts/config_phenos.yml v4/$tissue/phenos/config.yml
    # (comment out paired_end: True in config.yml for NAcc)
done

mkdir -p v4/NAcc/phenos/intermediate/{alt_TSS_polyA,expression,splicing,stability}
for tissue in NAcc1 NAcc2 NAcc3 NAcc4; do
    while read id; do
        for group in grp_1.downstream grp_1.upstream grp_2.downstream grp_2.upstream; do
            mkdir -p v4/NAcc/phenos/intermediate/alt_TSS_polyA/${group}
            ln -s ../../../../../${tissue}/phenos/intermediate/alt_TSS_polyA/${group}/${id} v4/NAcc/phenos/intermediate/alt_TSS_polyA/${group}/
        done
        ln -s ../../../../${tissue}/phenos/intermediate/expression/${id} v4/NAcc/phenos/intermediate/expression/
        ln -s ../../../../${tissue}/phenos/intermediate/splicing/${id}.junc v4/NAcc/phenos/intermediate/splicing/
        for type in exonic intronic; do
            ln -s ../../../../${tissue}/phenos/intermediate/stability/${id}.${type}.counts.txt v4/NAcc/phenos/intermediate/stability/
        done
    done < v4/${tissue}/rat_ids.txt
done

mkdir -p v4/PL/phenos/intermediate/{alt_TSS_polyA,expression,splicing,stability}
for tissue in PL1 PL2 PL3; do
    while read id; do
        for group in grp_1.downstream grp_1.upstream grp_2.downstream grp_2.upstream; do
            mkdir -p v4/PL/phenos/intermediate/alt_TSS_polyA/${group}
            ln -s ../../../../../${tissue}/phenos/intermediate/alt_TSS_polyA/${group}/${id} v4/PL/phenos/intermediate/alt_TSS_polyA/${group}/
        done
        ln -s ../../../../${tissue}/phenos/intermediate/expression/${id} v4/PL/phenos/intermediate/expression/
        ln -s ../../../../${tissue}/phenos/intermediate/splicing/${id}.junc v4/PL/phenos/intermediate/splicing/
        for type in exonic intronic; do
            ln -s ../../../../${tissue}/phenos/intermediate/stability/${id}.${type}.counts.txt v4/PL/phenos/intermediate/stability/
        done
    done < v4/${tissue}/rat_ids.txt
done

mkdir -p v4/pVTA/phenos/intermediate/{alt_TSS_polyA,expression,splicing,stability}
for tissue in pVTA1 pVTA2; do
    while read id; do
        for group in grp_1.downstream grp_1.upstream grp_2.downstream grp_2.upstream; do
            mkdir -p v4/pVTA/phenos/intermediate/alt_TSS_polyA/${group}
            ln -s ../../../../../${tissue}/phenos/intermediate/alt_TSS_polyA/${group}/${id} v4/pVTA/phenos/intermediate/alt_TSS_polyA/${group}/
        done
        ln -s ../../../../${tissue}/phenos/intermediate/expression/${id} v4/pVTA/phenos/intermediate/expression/
        ln -s ../../../../${tissue}/phenos/intermediate/splicing/${id}.junc v4/pVTA/phenos/intermediate/splicing/
        for type in exonic intronic; do
            ln -s ../../../../${tissue}/phenos/intermediate/stability/${id}.${type}.counts.txt v4/pVTA/phenos/intermediate/stability/
        done
    done < v4/${tissue}/rat_ids.txt
done
