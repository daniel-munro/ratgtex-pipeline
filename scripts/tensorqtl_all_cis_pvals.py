import argparse
from fastparquet import ParquetFile
import pandas as pd

p = argparse.ArgumentParser(description="Compile p-values for all tested cis-window variants.")
p.add_argument("nom_prefix", help="Prefix to nominal parquet files (e.g. {prefix}.cis_qtl_pairs.chr1.parquet)")
p.add_argument("outfile", help="Output file name (*.tsv.gz)")
args = p.parse_args()

df = []
for chrom in range(1, 21):
    f = f"{args.nom_prefix}.cis_qtl_pairs.chr{chrom}.parquet"
    d = ParquetFile(f).to_pandas()
    d = d[["phenotype_id", "variant_id", "pval_nominal"]]
    df.append(d)
df = pd.concat(df, ignore_index=True)
df.to_csv(args.outfile, sep="\t", index=False, float_format="%g")
