import argparse
from gtfparse import read_gtf
import numpy as np
import os
import pandas as pd
import rnaseqnorm

parser = argparse.ArgumentParser(description="Combine RSEM output files into BED files with log2 counts, TPM, and IQN units")
parser.add_argument("rsem_dir", help="Directory with RSEM output files to be combined. Each sample is named from its file prefix.")
parser.add_argument("annotations", help="GTF annotation file.")
parser.add_argument("out_prefix", help="Prefix (including path) of output files, which will be appended with '.log2.bed', '.tpm.bed', and '.iqn.bed'")
args = parser.parse_args()

files = os.listdir(args.rsem_dir)
ids = [x.strip(".genes.results.gz") for x in files]

for i, file in enumerate(files):
    d = pd.read_csv(os.path.join(args.rsem_dir, file), sep="\t", index_col="gene_id")
    if i == 0:
        log2 = pd.DataFrame(index=d.index)
        tpm = pd.DataFrame(index=d.index)
    else:
        assert d.index.equals(log2.index)
    log2[ids[i]] = np.log2(d["expected_count"] + 1)
    tpm[ids[i]] = d["TPM"]

iqn = rnaseqnorm.normalize_quantiles(tpm)
iqn = rnaseqnorm.inverse_normal_transform(iqn)

anno = read_gtf(args.annotations)
anno = anno.loc[anno["feature"] == "gene", :]
anno["tss"] = np.where(anno["strand"] == "+", anno["start"], anno["end"])
anno["start"] = anno["tss"] - 1  # BED coordinates are 0-based
anno["end"] = anno["tss"]
# anno = anno.loc[anno["seqname"].isin([str(i) for i in range(1, 21)]), :]
# anno["#chr"] = anno["seqname"].astype(int)  # for sorting
anno["#chr"] = anno["seqname"]
anno = anno.sort_values(["#chr", "start"])
anno = anno[["#chr", "start", "end", "gene_id"]]

# Save filtered IQN file to avoid tensorQTL error on phenotypes with 1 nonzero value:
iqnfilt = iqn.copy()
iqnfilt = iqnfilt.loc[np.sum(iqnfilt > 0, axis=1) > 1, :]

log2 = anno.merge(log2.reset_index(), on="gene_id", how="inner")
tpm = anno.merge(tpm.reset_index(), on="gene_id", how="inner")
iqn = anno.merge(iqn.reset_index(), on="gene_id", how="inner")
iqnfilt = anno.merge(iqnfilt.reset_index(), on="gene_id", how="inner")

iqnfilt = iqnfilt.loc[iqnfilt["#chr"].isin([str(i) for i in range(1, 21)]), :]
iqnfilt["#chr"] = iqnfilt["#chr"].astype(int)  # for sorting
iqnfilt = iqnfilt.sort_values(["#chr", "start"])

log2.to_csv(f"{args.out_prefix}.log2.bed", sep="\t", index=False, float_format="%g")
tpm.to_csv(f"{args.out_prefix}.tpm.bed", sep="\t", index=False, float_format="%g")
iqn.to_csv(f"{args.out_prefix}.iqn.bed", sep="\t", index=False, float_format="%g")
iqnfilt.to_csv(f"{args.out_prefix}.iqn.filtered.bed", sep="\t", index=False, float_format="%g")
