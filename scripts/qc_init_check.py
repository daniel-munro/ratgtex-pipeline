"""Scan raw data and setup to report stats and issues"""

import argparse
from pathlib import Path
import pysam
import yaml

def load_fastq_lists(map_file: Path, fastq_dir: Path, paired_end=None) -> dict:
    """Load the FASTQ file paths for each sample.
    
    Returns a dictionary where each value is a tuple (paths, is_paired).
    paths is a list of all single-end or paired-end FASTQ files for the sample.
    is_paired is a boolean indicating if the sample is paired-end.
    """
    paths = {}
    with open(map_file, "r") as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            rat_id = fields[-1]
            is_paired = len(fields) == 3
            
            if rat_id not in paths:
                paths[rat_id] = ([], is_paired)
            elif paths[rat_id][1] != is_paired:
                raise ValueError(f"Sample {rat_id} has mixed paired-end and single-end reads")
            
            paths[rat_id][0].append(fastq_dir / fields[0])
            if is_paired:
                paths[rat_id][0].append(fastq_dir / fields[1])
    # Validate against config if specified
    if paired_end is not None:
        for rat_id, (_, is_paired) in paths.items():
            if is_paired != paired_end:
                raise ValueError(f"Sample {rat_id} has {'paired' if is_paired else 'single'}-end reads, but config specifies {'paired' if paired_end else 'single'}-end")
    
    return paths

parser = argparse.ArgumentParser(description="Scan raw data and setup to report stats and issues")
parser.add_argument("version", type=str, help="RatGTEx version, e.g. v1")
parser.add_argument("tissue", type=str, help="Tissue name")
args = parser.parse_args()

stats = {}

with open(f"{args.version}/{args.tissue}/rat_ids.txt", "r") as f:
    ids = f.read().splitlines()

config = yaml.safe_load(open("config.yaml"))["tissues"][args.tissue]
read_length = config["read_length"]
fastq_path = Path(config["fastq_path"])
paired = config.get("paired_end", None)
geno_dataset = config["geno_dataset"]

print(f"INFO: Paired-end configuration: {paired}")
print(f"INFO: Read length is set to {read_length} in config.yaml. Ensure this is correct.")
print(f"INFO: Using geno/{geno_dataset}.vcf.gz")

map_file = f"{args.version}/{args.tissue}/fastq_map.txt"
fastqs = load_fastq_lists(map_file, fastq_path, paired)

for rat_id in ids:
    if not rat_id in fastqs:
        raise ValueError(f"{rat_id} not found in FASTQ map")

print("PASS: FASTQ map contains all individuals to be included")

n_fastq_ids = len(fastqs.keys())
print(f"INFO: {len(ids)} out of {n_fastq_ids} individuals in FASTQ map will be used")

for rat_id in ids:
    paths, paired = fastqs[rat_id]
    for fastq in paths:
        if not fastq.exists():
            raise FileNotFoundError(f"{fastq} not found")
print("PASS: All necessary FASTQ files exist")

vcf = pysam.VariantFile(f"geno/{geno_dataset}.vcf.gz")
samples = list(vcf.header.samples)
missing = [id for id in ids if not id in samples]
if len(missing) > 0:
    raise ValueError(f"{len(missing)} individuals missing from VCF: {missing}")
print("PASS: VCF contains all individuals to be included")
