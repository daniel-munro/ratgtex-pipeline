"""Compare ASE read counts to genotypes to check for mismatches"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import pysam


def geno_frac_alt(gt: tuple) -> float:
    """Convert genotype to fraction alt allele (0, 0.5, or 1)"""
    return None if gt[0] is None or gt[1] is None else sum(gt) / 2


parser = argparse.ArgumentParser(
    description="Compare ASE read counts to genotypes to check for mismatches"
)
parser.add_argument("vcf", type=Path, help="VCF file used as ASEReadCounter input")
parser.add_argument(
    "count_dir",
    type=Path,
    help="Directory containing ASEReadCounter outputs named {sample_id}.readcounts.txt",
)
parser.add_argument("samples", type=Path, help="File containing sample IDs, one per line, to check")
parser.add_argument(
    "out_matrix", type=Path, help="Output file for matrix of RNA/genotype similarities"
)
parser.add_argument("out_summary", type=Path, help="Output file for summary table")
args = parser.parse_args()

with args.samples.open() as f:
    samples = f.read().splitlines()
vcf = pysam.VariantFile(args.vcf)
geno_samples = list(vcf.header.samples)
geno = {sample: {} for sample in geno_samples}
for rec in vcf.fetch():
    for sample in geno_samples:
        geno[sample][rec.id] = geno_frac_alt(rec.samples[sample]["GT"])

# Make table of similarities for RNA samples vs. genotype rats
simil = pd.DataFrame(index=samples, columns=geno_samples, dtype=float)
# Track number of variants used for each comparison
variant_counts = pd.DataFrame(index=samples, columns=geno_samples, dtype=int)

for sample in samples:
    counts = pd.read_csv(
        args.count_dir / f"{sample}.readcounts.txt", sep="\t", index_col=2
    )
    counts["frac_alt"] = counts["altCount"] / (counts["refCount"] + counts["altCount"])
    for geno_sam in geno_samples:
        snps = [snp for snp in counts.index if snp in geno[geno_sam] and geno[geno_sam][snp] is not None]
        diff = counts.loc[snps, "frac_alt"] - [geno[geno_sam][snp] for snp in snps]
        sim = np.mean(1 - np.abs(diff))
        simil.loc[sample, geno_sam] = sim
        variant_counts.loc[sample, geno_sam] = len(snps)

simil.to_csv(args.out_matrix, sep="\t", index_label="RNA_ID", float_format="%g")

# Make table of top similarities per RNA sample for quick checking
# Also include 2nd highest to show how big the gap is
top = pd.DataFrame(index=samples)
top["top_simil_ID"] = [
    simil.columns[np.argmax(simil.loc[sample, :])] for sample in top.index
]
top["top_simil"] = simil.max(axis=1)
top["top_variant_count"] = [
    variant_counts.loc[sample, top.loc[sample, "top_simil_ID"]] for sample in top.index
]
top["is_correct"] = top.index == top["top_simil_ID"]
top["next_simil_ID"] = [simil.loc[sample].nlargest(2).index[-1] for sample in top.index]
top["next_simil"] = [simil.loc[sample].nlargest(2).iloc[-1] for sample in top.index]
top.to_csv(args.out_summary, sep="\t", index_label="RNA_ID", float_format="%g")
