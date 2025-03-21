"""Scan raw data and setup to report stats and issues"""

import argparse
import pandas as pd
from pathlib import Path
import pysam
import yaml

parser = argparse.ArgumentParser(description="Scan raw data and setup to report stats and issues")
parser.add_argument("version", type=str, help="RatGTEx version, e.g. v1")
parser.add_argument("tissue", type=str, help="Tissue name")
args = parser.parse_args()

stats = {}

with open(f"{args.version}/{args.tissue}/rat_ids.txt", "r") as f:
    ids = f.read().splitlines()

config = yaml.safe_load(open("config.yaml"))[args.tissue]
read_length = config["read_length"]
fastq_path = Path(config["fastq_path"])
paired = bool(config["paired_end"])
geno_dataset = config["geno_dataset"]

print(f"INFO: Paired-end sequencing: {paired}")
print(f"INFO: Read length is set to {read_length} in config.yaml. Ensure this is correct.")
print(f"INFO: Using geno/{geno_dataset}.vcf.gz")

names = ["fastq1", "fastq2", "rat_id"] if paired else ["fastq", "rat_id"]
fastqs = pd.read_csv(f"{args.version}/{args.tissue}/fastq_map.txt", sep="\t", names=names, dtype=str)

for rat_id in ids:
    if not rat_id in fastqs["rat_id"].values:
        raise ValueError(f"{rat_id} not found in FASTQ map")

print("PASS: FASTQ map contains all individuals to be included")

n_fastq_ids = len(fastqs["rat_id"].unique())
print(f"INFO: {len(ids)} out of {n_fastq_ids} individuals in FASTQ map will be used")

fastqs = fastqs.loc[fastqs["rat_id"].isin(ids)]
for i in range(fastqs.shape[0]):
    if paired:
        fastq1 = fastq_path / fastqs["fastq1"].iloc[i]
        if not fastq1.exists():
            raise FileNotFoundError(f"{fastq1} not found")
        fastq2 = fastq_path / fastqs["fastq2"].iloc[i]
        if not fastq2.exists():
            raise FileNotFoundError(f"{fastq2} not found")
    else:
        fastq = fastq_path / fastqs["fastq"].iloc[i]
        if not fastq.exists():
            raise FileNotFoundError(f"{fastq} not found")
print("PASS: All necessary FASTQ files exist")

vcf = pysam.VariantFile(f"geno/{geno_dataset}.vcf.gz")
samples = list(vcf.header.samples)
missing = [id for id in ids if not id in samples]
if len(missing) > 0:
    raise ValueError(f"{len(missing)} individuals missing from VCF: {missing}")
print("PASS: VCF contains all individuals to be included")
