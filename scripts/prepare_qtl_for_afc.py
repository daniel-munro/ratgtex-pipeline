# Extract from tensorQTL output file the phenotype_id and
# variant_id column. Then extract chromosome and location
# from variant IDs and include those too. You can include multiple file names,
# and the union of gene/SNP pairs will be output. Useful to get aFC for top associations +
# independent eQTLs.

import sys
import pandas as pd

d = []
for infile in sys.argv[1:]:
    df = pd.read_csv(infile, sep="\t")
    df = df[["phenotype_id", "variant_id"]]
    d.append(df)
d = pd.concat(d)
d = d.drop_duplicates()
print("pid\tsid\tsid_chr\tsid_pos")
for row in d.to_dict(orient="records"):
    chrom, pos = row["variant_id"].split(":")
    print(f"{row['phenotype_id']}\t{row['variant_id']}\t{chrom}\t{pos}")
