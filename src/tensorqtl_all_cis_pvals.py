import argparse
from fastparquet import ParquetFile
from glob import glob
import pandas as pd

p = argparse.ArgumentParser(description="Compile p-values for all tested cis-window variants.")
p.add_argument("parquets", help="Directory containing parquet files with all tested gene-variant pairs")
p.add_argument("outfile", help="Output file name (*.txt.gz)")
args = p.parse_args()

df = []
for f in glob(f"{args.parquets}/*.parquet"):
    d = ParquetFile(f).to_pandas()
    d = d[["phenotype_id", "variant_id", "pval_nominal"]]
    df.append(d)
df = pd.concat(df)
df.to_csv(args.outfile, sep="\t", index=False, float_format="%g")
