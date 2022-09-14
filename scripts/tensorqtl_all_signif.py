import pandas as pd
import numpy as np
import argparse
from fastparquet import ParquetFile

parser = argparse.ArgumentParser(
    description="Get per-gene pval thresholds from tensorQTL cis (permutations) output, then use it to filter cis_nominal output to get all significant cis-eQTLs."
)
parser.add_argument("perm_file", help="tensorQTL perm file (*.cis_qtl.txt.gz)")
parser.add_argument("nom_prefix", help="prefix to nominal parquet files (e.g. {prefix}.cis_qtl_pairs.1.parquet)")
parser.add_argument("output", help="output file (TSV)")
args = parser.parse_args()

perm = pd.read_csv(args.perm_file, sep="\t")
# cutoff = {}
# for i in range(perm.shape[0]):
#     cutoff[perm.phenotype_id[i]] = perm.pval_nominal_threshold[i]
perm = perm[["phenotype_id", "pval_nominal_threshold"]]

sig = []
# Read each nominal file and save significant pairs.
nom = []
for i in range(1, 21):
    d = ParquetFile(f"{args.nom_prefix}.cis_qtl_pairs.{i}.parquet").to_pandas()
    d = d.merge(perm, how="inner", on="phenotype_id")
    # print(d.loc[d.pval_nominal_threshold.isnull(), :])
    assert np.sum(d.pval_nominal_threshold.isnull()) == 0
    d = d.loc[d.pval_nominal < d.pval_nominal_threshold, :]
    sig.append(d)
sig = pd.concat(sig, ignore_index=True)
# sig = sig.rename(columns={"phenotype_id": "gene_id"})
# sig = sig.drop(columns=["tss_distance"])

sig.to_csv(args.output, sep="\t", index=False, float_format="%g")
