"""Map requested biological rats to PLINK samples without guessing among replicates."""
from collections import defaultdict
from pathlib import Path

import pandas as pd


def metadata_sample_ids(path):
    """Construct PSAM join keys from declared metadata fields, without parsing IDs."""
    metadata = pd.read_csv(path, dtype=str, keep_default_na=False,
                           usecols=['rfid', 'library_name', 'barcode'])
    if metadata[['rfid', 'library_name']].eq('').any().any():
        raise ValueError('Genotype metadata has missing RFID or library_name')
    plain = metadata.library_name + '_' + metadata.rfid
    candidates = pd.concat([
        pd.DataFrame({'IID': plain, 'rfid': metadata.rfid}),
        pd.DataFrame({'IID': plain + '_' + metadata.barcode, 'rfid': metadata.rfid})[
            metadata.barcode.ne('')],
    ], ignore_index=True).drop_duplicates()
    if candidates.IID.duplicated().any():
        conflicts = candidates.loc[candidates.IID.duplicated(False), 'IID'].unique()
        raise ValueError('Conflicting metadata RFIDs for genotype ID(s): ' + ', '.join(conflicts))
    return candidates


def prepare_samples(sample_files, datasets, psams, metadata, keep, rename, coverage, counts):
    frames = [pd.read_csv(path, sep=r'\s+', dtype=str).rename(
        columns={'#IID': 'IID', '#FID': 'FID'}) for path in psams]
    first = frames[0]
    if first.IID.duplicated().any():
        raise ValueError('Duplicate genotype IIDs')
    for path, frame in zip(psams[1:], frames[1:]):
        if not first.equals(frame):
            raise ValueError(f'Chromosome PSAM differs from {psams[0]}: {path}')
    annotated = first.merge(metadata_sample_ids(metadata), on='IID', how='left', validate='one_to_one')
    if annotated.rfid.isna().any():
        missing = annotated.loc[annotated.rfid.isna(), 'IID']
        raise ValueError('PSAM sample(s) absent from genotype metadata: ' + ', '.join(missing))
    mapping = defaultdict(list)
    for row in annotated.to_dict('records'):
        mapping[row['rfid']].append((row.get('FID', '0'), row['IID']))
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
                    snakemake.input.psams, snakemake.input.metadata, snakemake.output.keep, snakemake.output.rename,
                    snakemake.log.coverage, snakemake.log.counts)
