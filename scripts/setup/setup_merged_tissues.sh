set -e

## fastq_map.txt are created just to provide for download

mkdir -p v3/NAcc/{rsem_out,splice/junc}
cat v3/NAcc{1,2,3}/rat_ids.txt | sort | uniq > v3/NAcc/rat_ids.txt

# Prepend sub-tissue ID to each FASTQ path
awk '{for(i=1;i<NF;i++) printf "NAcc1/%s\t", $i; print $NF}' v3/NAcc1/fastq_map.txt > v3/NAcc/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "NAcc2/%s\t", $i; print $NF}' v3/NAcc2/fastq_map.txt >> v3/NAcc/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "NAcc3/%s\t", $i; print $NF}' v3/NAcc3/fastq_map.txt >> v3/NAcc/fastq_map.txt

for tissue in NAcc1 NAcc2 NAcc3; do
    while read id; do
        ln -s ../../${tissue}/rsem_out/${id}.genes.results.gz v3/NAcc/rsem_out/
        ln -s ../../../${tissue}/splice/junc/${id}.junc.gz v3/NAcc/splice/junc/
    done < v3/${tissue}/rat_ids.txt
done

mkdir -p v3/PL/{rsem_out,splice/junc}
cat v3/PL{1,2,3}/rat_ids.txt | sort | uniq > v3/PL/rat_ids.txt

# Prepend sub-tissue ID to each FASTQ path
awk '{for(i=1;i<NF;i++) printf "PL1/%s\t", $i; print $NF}' v3/PL1/fastq_map.txt > v3/PL/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "PL2/%s\t", $i; print $NF}' v3/PL2/fastq_map.txt >> v3/PL/fastq_map.txt
awk '{for(i=1;i<NF;i++) printf "PL3/%s\t", $i; print $NF}' v3/PL3/fastq_map.txt >> v3/PL/fastq_map.txt

for tissue in PL1 PL2 PL3; do
    while read id; do
        ln -s ../../${tissue}/rsem_out/${id}.genes.results.gz v3/PL/rsem_out/
        ln -s ../../../${tissue}/splice/junc/${id}.junc.gz v3/PL/splice/junc/
    done < v3/${tissue}/rat_ids.txt
done
