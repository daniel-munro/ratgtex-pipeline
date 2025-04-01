"""Predict sex from chrY gene expression and compare to metadata labels"""

import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Predict sex from chrY gene expression and compare to metadata labels"
)
parser.add_argument("expr", type=str, help="Expression BED file with TPM values")
parser.add_argument("meta", type=str, help="Metadata, currently expects a CSV table that includes columns 'rfid' and 'sex'")
parser.add_argument("out", type=str, help="Output file for report summarizing sex concordance")
args = parser.parse_args()

# Cutoff for total chrY gene TPM, empirically found to predict female/male
TPM_CUTOFF = 4
# Cap individual TPM values before summing to control for outliers
MAX_TPM = 5

expr = pd.read_csv(args.expr, sep="\t", dtype={"#chr": str})
expr = expr.loc[expr["#chr"] == "chrY"]
expr = expr.drop(columns=["#chr", "start", "end", "gene_id"])
expr = expr.clip(upper=MAX_TPM)
expr = expr.sum(axis=0) # Sum TPM for chrY genes per sample
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
expr = expr.join(meta, how="inner")

out = open(args.out, "w")
out.write(f"Sex labels present for {len(expr.index)} of {n_samples} samples.\n")
if len(expr.index) > 0:
    # Count matches between pred and label columns:
    n_matches = (expr["pred"] == expr["label"]).sum()
    out.write(f"{n_matches} matches between expression-predicted sex and metadata labels.\n")
    mismatches = expr.loc[expr["pred"] != expr["label"]]
    if len(mismatches.index) > 0:
        out.write(f"Samples with mismatch (chrY total TPM cutoff = {TPM_CUTOFF}):\n")
        out.write(mismatches.to_string() + "\n")
    out.write("\nAll results:\n")
    out.write(expr.to_string() + "\n")
out.close()
