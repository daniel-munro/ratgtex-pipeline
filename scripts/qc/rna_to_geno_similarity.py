"""Compare ASE read counts to genotypes to check for mismatches."""

import argparse
import logging
from pathlib import Path
from time import monotonic
import numpy as np
import pandas as pd
import pysam


def geno_frac_alt(gt: tuple) -> float:
    """Convert diploid biallelic calls to alt fraction; exclude other calls."""
    if len(gt) == 2 and gt[0] in {0, 1} and gt[1] in {0, 1}:
        return sum(gt) / 2
    return None


parser = argparse.ArgumentParser(
    description="Compare ASE read counts to genotypes to check for mismatches"
)
parser.add_argument("--vcf", type=Path, required=True, help="VCF containing candidate genotypes")
parser.add_argument(
    "--count-dir", required=True,
    type=Path,
    help="Directory containing ASEReadCounter outputs named {sample_id}.readcounts.txt",
)
parser.add_argument("--samples", type=Path, required=True, help="File containing RNA sample IDs, one per line, to check")
parser.add_argument(
    "--out-matrix", type=Path, help="Optional output file for matrix of RNA/genotype similarities"
)
parser.add_argument("--out-summary", type=Path, required=True, help="Output file for summary table")
parser.add_argument("--snps", type=Path, help="Variant IDs, one per line, from which --n-snps samples")
parser.add_argument("--n-snps", type=int, help="Randomly select N variant IDs using seed 20230312")
args = parser.parse_args()
if args.n_snps is not None and (args.snps is None or args.n_snps <= 0):
    parser.error("--n-snps requires --snps and a positive SNP count")
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
started = monotonic()

with args.samples.open() as f:
    samples = f.read().splitlines()
snp_set = None
if args.n_snps is not None:
    with args.snps.open() as f:
        snps = f.read().splitlines()
    rng = np.random.default_rng(seed=20230312)
    snp_set = set(rng.choice(snps, size=args.n_snps, replace=False))
    logging.info("Selected %d variant IDs using seed 20230312", len(snp_set))
vcf = pysam.VariantFile(args.vcf)
geno_samples = list(vcf.header.samples)
logging.info("Loading VCF genotypes for %d candidates", len(geno_samples))
variant_rows = {}
geno_rows = []
for scanned, rec in enumerate(vcf.fetch(), start=1):
    if scanned % 100000 == 0:
        logging.info("Scanned %d VCF records; loaded %d variants (%.1f s elapsed)",
                     scanned, len(geno_rows), monotonic() - started)
    if snp_set is not None and rec.id not in snp_set:
        continue
    row = len(geno_rows)
    # Retain the last record for repeated IDs, as the original dictionaries did.
    variant_rows[rec.id] = row
    geno_rows.append(np.fromiter(
        (geno_frac_alt(call["GT"]) for call in rec.samples.values()),
        dtype=np.float32, count=len(geno_samples),
    ))
    if (row + 1) % 5000 == 0:
        logging.info("Loaded %d variants (%.1f s elapsed)", row + 1, monotonic() - started)
geno = np.stack(geno_rows)
del geno_rows
vcf.close()
logging.info("Loaded %d variants; comparing %d RNA samples", len(geno), len(samples))

# Make table of similarities for RNA samples vs. genotype rats
similarities = np.empty((len(samples), len(geno_samples)), dtype=float)
# Track number of variants used for each comparison
comparison_counts = np.empty(similarities.shape, dtype=np.int64)

for sample_row, sample in enumerate(samples):
    counts = pd.read_csv(
        args.count_dir / f"{sample}.readcounts.txt", sep="\t", index_col=2
    )
    counts["frac_alt"] = counts["altCount"] / (counts["refCount"] + counts["altCount"])
    # Align SNP IDs once per RNA sample, then compare candidate blocks in NumPy.
    matched = counts.index.isin(variant_rows)
    rows = np.array([variant_rows[snp] for snp in counts.index[matched]], dtype=np.int64)
    fractions = counts.loc[matched, "frac_alt"].to_numpy()
    for start in range(0, len(geno_samples), 128):
        stop = min(start + 128, len(geno_samples))
        candidate_genotypes = geno[rows, start:stop]
        comparison_counts[sample_row, start:stop] = np.sum(
            ~np.isnan(candidate_genotypes), axis=0
        )
        scores = 1 - np.abs(fractions[:, None] - candidate_genotypes)
        # pandas mean skips NaN RNA fractions too, but variant counts include them.
        valid = np.sum(~np.isnan(scores), axis=0)
        similarities[sample_row, start:stop] = np.divide(
            np.nansum(scores, axis=0), valid,
            out=np.full(stop - start, np.nan), where=valid > 0,
        )
    if (sample_row + 1) % 10 == 0 or sample_row + 1 == len(samples):
        logging.info("Compared %d/%d RNA samples (%.1f s elapsed)",
                     sample_row + 1, len(samples), monotonic() - started)

simil = pd.DataFrame(similarities, index=samples, columns=geno_samples)
variant_counts = pd.DataFrame(comparison_counts, index=samples, columns=geno_samples)

if args.out_matrix is not None:
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
logging.info("Wrote summary%s (%.1f s elapsed)",
             " and similarity matrix" if args.out_matrix is not None else "",
             monotonic() - started)
