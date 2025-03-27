set -e

mkdir -p v3/NAcc/{rsem_out,splice/junc}
cat v3/NAcc{1,2,3}/rat_ids.txt | sort | uniq > v3/NAcc/rat_ids.txt
for tissue in NAcc1 NAcc2 NAcc3; do
    while read id; do
        ln -s ../../${tissue}/rsem_out/${id}.genes.results.gz v3/NAcc/rsem_out/
        ln -s ../../../${tissue}/splice/junc/${id}.junc.gz v3/NAcc/splice/junc/
    done < v3/${tissue}/rat_ids.txt
done

mkdir -p v3/PL/{rsem_out,splice/junc}
cat v3/PL{1,2,3}/rat_ids.txt | sort | uniq > v3/PL/rat_ids.txt
for tissue in PL1 PL2 PL3; do
    while read id; do
        ln -s ../../${tissue}/rsem_out/${id}.genes.results.gz v3/PL/rsem_out/
        ln -s ../../../${tissue}/splice/junc/${id}.junc.gz v3/PL/splice/junc/
    done < v3/${tissue}/rat_ids.txt
done
