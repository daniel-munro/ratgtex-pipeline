"""Map requested biological rats to PLINK samples without guessing among replicates."""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

PATTERNS = [
    r'^ddGBSR?_[^_]+_[^_]+_(.+)$',
    r'^(?:Riptide|Rattaca)[^_]+_(.+)_[^_]+$',
    r'^LSWrnaseq_(.+)_LSWNNNNN$',
]


def prepare_samples(sample_files, datasets, psams, keep, rename, coverage, counts):
    frames = [pd.read_csv(path, sep=r'\s+', dtype=str).rename(
        columns={'#IID': 'IID', '#FID': 'FID'}) for path in psams]
    first = frames[0]
    if first.IID.duplicated().any():
        raise ValueError('Duplicate genotype IIDs')
    for path, frame in zip(psams[1:], frames[1:]):
        if not first.equals(frame):
            raise ValueError(f'Chromosome PSAM differs from {psams[0]}: {path}')
    mapping = defaultdict(list)
    for row in first.to_dict('records'):
        for pattern in PATTERNS:
            match = re.fullmatch(pattern, row['IID'])
            if match:
                mapping[match[1]].append((row.get('FID', '0'), row['IID']))
                break
    records = []
    for dataset, path in zip(datasets, sample_files):
        rfids = Path(path).read_text().splitlines()
        if not rfids or any(not r or r != r.strip() for r in rfids) or len(set(rfids)) != len(rfids):
            raise ValueError(f'Empty, duplicate, or malformed requested RFIDs in {path}')
        for rfid in rfids:
            matches = mapping[rfid]
            status = 'present' if len(matches) == 1 else ('missing' if not matches else 'ambiguous')
            records.append((dataset, rfid, status, ';'.join(x[1] for x in matches)))
    table = pd.DataFrame(records, columns=['dataset', 'rfid', 'status', 'genotype_ids'])
    table.to_csv(coverage, sep='\t', index=False)
    totals = table.groupby(['dataset', 'status']).size().unstack(fill_value=0).reindex(
        columns=['present', 'missing', 'ambiguous'], fill_value=0)
    totals.insert(0, 'requested', totals.sum(axis=1))
    totals.to_csv(counts, sep='\t')
    if (table.status != 'present').any():
        raise ValueError(f'Missing or ambiguous genotype identities; inspect {coverage} and {counts}')
    rfids = sorted(set(table.rfid))
    with open(keep, 'w') as out_keep, open(rename, 'w') as out_rename:
        out_keep.write('#FID\tIID\n')
        out_rename.write('#OLD_FID\tOLD_IID\tNEW_FID\tNEW_IID\n')
        for rfid in rfids:
            fid, iid = mapping[rfid][0]
            out_keep.write(f'{fid}\t{iid}\n')
            out_rename.write(f'{fid}\t{iid}\t0\t{rfid}\n')


if __name__ == '__main__':
    prepare_samples(snakemake.input.samples, snakemake.params.datasets,
                    snakemake.input.psams, snakemake.output.keep, snakemake.output.rename,
                    snakemake.log.coverage, snakemake.log.counts)
