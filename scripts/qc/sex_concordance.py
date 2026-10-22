"""Predict sex from chrY gene expression and compare to metadata labels"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Predict sex from chrY gene expression and compare to metadata labels"
)
parser.add_argument("expr", type=str, help="Expression BED file with TPM values")
parser.add_argument("meta", type=str, help="Metadata, currently expects a CSV table that includes columns 'rfid' and 'sex'")
parser.add_argument("out", type=str, help="Output file for report summarizing sex concordance")
parser.add_argument("--tpm-cutoff", type=float, default=30,
                    help="Male if capped chrY TPM sum >= cutoff (default: 30; calibrated on v5 Eye)")
parser.add_argument("--max-tpm", type=float, default=5,
                    help="Cap each gene TPM before summing (default: 5)")
args = parser.parse_args()

# Provisional v5 cutoff: Eye females 6.66-26.04, males 33.72-58.17.
# Reassess separation in other tissues rather than fitting to their labels.
TPM_CUTOFF = args.tpm_cutoff
# Cap individual TPM values before summing to control for outliers
MAX_TPM = args.max_tpm

expr = pd.read_csv(args.expr, sep="\t", dtype={"#chr": str})
expr = expr.loc[expr["#chr"] == "chrY"]
n_y_genes = len(expr)
if n_y_genes == 0:
    raise ValueError("Expression table contains no chrY genes")
expr = expr.drop(columns=["#chr", "start", "end", "phenotype_id"])
if expr.isna().any().any() or (expr < 0).any().any():
    raise ValueError("chrY TPM values must be nonmissing and nonnegative")
expr = expr.clip(upper=MAX_TPM)
expr = expr.sum(axis=0, skipna=False) # Sum TPM for chrY genes per sample
# Convert to DataFrame with sample names as index
expr = pd.DataFrame(expr, columns=["TPM"])
# Predict sex as F if TPM < cutoff, M if TPM >= cutoff:
expr["pred"] = expr["TPM"].apply(lambda x: "F" if x < TPM_CUTOFF else "M")

meta = pd.read_csv(args.meta, sep=",", dtype=str)
# Keep only the last entry per sample
meta = meta.drop_duplicates(subset=['rfid'], keep='last')
meta = meta.set_index("rfid")
meta = meta[["sex"]].rename(columns={"sex": "label"})
n_samples = len(expr.index)
expr = expr.join(meta, how="left")
labeled = expr.loc[expr["label"].isin(["F", "M"])]

out = open(args.out, "w")
out.write(f"chrY genes: {n_y_genes}; per-gene TPM cap: {MAX_TPM}; male cutoff: {TPM_CUTOFF}.\n")
out.write(f"F/M sex labels present for {len(labeled)} of {n_samples} samples.\n")
out.write("\nCapped chrY TPM sum by metadata sex:\n")
out.write(labeled.groupby("label")["TPM"].agg(["count", "min", "median", "max"]).to_string() + "\n\n")
if len(labeled.index) > 0:
    # Count matches between pred and label columns:
    n_matches = (labeled["pred"] == labeled["label"]).sum()
    out.write(f"{n_matches} matches between expression-predicted sex and metadata labels.\n")
    mismatches = labeled.loc[labeled["pred"] != labeled["label"]]
    if len(mismatches.index) > 0:
        out.write(f"Samples with mismatch (chrY total TPM cutoff = {TPM_CUTOFF}):\n")
        out.write(mismatches.to_string() + "\n")
out.write("\nAll results (TPM = capped chrY sum):\n")
out.write(expr.to_string() + "\n")
out.close()
