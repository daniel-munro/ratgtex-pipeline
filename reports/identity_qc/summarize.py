"""Collect identity mismatches and uninformative results from configured datasets."""

from pathlib import Path

import pandas as pd

records = []
for tissue, datasets in snakemake.config['tissues'].items():
    for dataset in datasets:
        path = Path(snakemake.config['version']) / dataset / 'qc/rna_to_geno_summary.tsv'
        frame = pd.read_csv(path, sep='\t', dtype={'RNA_ID': str, 'top_simil_ID': str,
                                                  'next_simil_ID': str})
        frame.insert(0, 'dataset', dataset)
        frame.insert(0, 'tissue', tissue)
        frame['similarity_gap'] = frame.top_simil - frame.next_simil
        reasons = []
        for row in frame.itertuples():
            issues = []
            if pd.isna(row.top_simil) or pd.isna(row.top_simil_ID):
                issues.append('no_similarity')
            elif row.RNA_ID != row.top_simil_ID:
                issues.append('identity_mismatch')
            if pd.isna(row.top_variant_count) or row.top_variant_count <= 0:
                issues.append('no_usable_variants')
            if pd.notna(row.similarity_gap) and row.similarity_gap <= 0:
                issues.append('ambiguous_top_match')
            reasons.append(';'.join(issues))
        frame['issue'] = reasons
        records.append(frame.loc[frame.issue != ''])
report = pd.concat(records, ignore_index=True)
report.to_csv(snakemake.output[0], sep='\t', index=False, float_format='%g')
print(f'{len(report)} identity-QC issues across {len(records)} datasets')
