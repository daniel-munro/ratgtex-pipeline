import argparse
from gtfparse import read_gtf
import numpy as np
import pandas as pd
from pathlib import Path
import rnaseqnorm

parser = argparse.ArgumentParser(description="Combine RSEM output files into BED files with log2 counts, TPM, and IQN units")
parser.add_argument("rsem_dir", type=Path, help="Directory with RSEM output files to be combined. Each sample is named from its file prefix.")
parser.add_argument("samples", type=Path, help="File containing sample IDs, one per line, for which to load RSEM outputs")
parser.add_argument("annotations", help="GTF annotation file.")
parser.add_argument("out_prefix", help="Prefix (including path) of output files, which will be appended with '.log2.bed', '.tpm.bed', and '.iqn.bed'")
args = parser.parse_args()

with args.samples.open() as f:
    ids = f.read().splitlines()

for i, sample in enumerate(ids):
    fname = args.rsem_dir / f"{sample}.genes.results.gz"
    d = pd.read_csv(fname, sep="\t", index_col="gene_id")
    if i == 0:
        log2 = pd.DataFrame(index=d.index)
        tpm = pd.DataFrame(index=d.index)
    else:
        assert d.index.equals(log2.index)
    log2[sample] = np.log2(d["expected_count"] + 1)
    tpm[sample] = d["TPM"]

iqn = rnaseqnorm.normalize_quantiles(tpm)
iqn = rnaseqnorm.inverse_normal_transform(iqn)

anno = read_gtf(args.annotations)
# Newer versions return a polars DF by default, but not all versions allow
# return type to be specified, so this handles older and newer versions:
if type(anno).__module__ == 'polars.dataframe.frame':
    anno = anno.to_pandas()
anno = anno.loc[anno["feature"] == "gene", :]
anno["tss"] = np.where(anno["strand"] == "+", anno["start"], anno["end"])
anno["start"] = anno["tss"] - 1  # BED coordinates are 0-based
anno["end"] = anno["tss"]
# anno = anno.loc[anno["seqname"].isin([str(i) for i in range(1, 21)]), :]
# anno["#chr"] = anno["seqname"].astype(int)  # for sorting
anno["#chr"] = anno["seqname"]
anno = anno.sort_values(["#chr", "start"])
anno = anno[["#chr", "start", "end", "gene_id"]]

# Save filtered IQN file to use for tensorQTL:
iqnfilt = iqn.copy()
# iqnfilt = iqnfilt.loc[np.sum(iqnfilt > 0, axis=1) > 1, :]

# Remove rows with more than 50% zeros, as recommended for tensorQTL:
mostly_zeros = tpm.apply(lambda x: np.mean(x == 0) > 0.5, axis=1)
assert tpm.index.equals(iqnfilt.index)
iqnfilt = iqnfilt[~mostly_zeros]
print(f"{args.out_prefix}: Removed {mostly_zeros.sum()} rows with more than 50% zeros")

# Remove rows with no variation:
no_var = iqnfilt.apply(lambda x: x.nunique() == 1, axis=1)
iqnfilt = iqnfilt[~no_var]
print(f"{args.out_prefix}: Removed {no_var.sum()} additional rows with no variation")

# Stop if the processed table is empty:
assert not iqnfilt.empty, f"After filtering, no phenotypes remain"

log2 = anno.merge(log2.reset_index(), on="gene_id", how="inner")
tpm = anno.merge(tpm.reset_index(), on="gene_id", how="inner")
iqn = anno.merge(iqn.reset_index(), on="gene_id", how="inner")
iqnfilt = anno.merge(iqnfilt.reset_index(), on="gene_id", how="inner")

# Filter to keep only chr1-chr20 and sort numerically by chromosome number
iqnfilt = iqnfilt.loc[iqnfilt["#chr"].isin([f"chr{i}" for i in range(1, 21)]), :]
# Extract chromosome number, convert to int for proper numeric sorting
iqnfilt["chr_num"] = iqnfilt["#chr"].str.extract(r"chr(\d+)").astype(int)
iqnfilt = iqnfilt.sort_values(["chr_num", "start"])
iqnfilt = iqnfilt.drop("chr_num", axis=1)

log2.to_csv(f"{args.out_prefix}.log2.bed", sep="\t", index=False, float_format="%g")
tpm.to_csv(f"{args.out_prefix}.tpm.bed", sep="\t", index=False, float_format="%g")
iqn.to_csv(f"{args.out_prefix}.iqn.bed", sep="\t", index=False, float_format="%g")
iqnfilt.to_csv(f"{args.out_prefix}.iqn.filtered.bed", sep="\t", index=False, float_format="%g")
