"""Join archival RNA-seq identities and prepare v5 inputs without exclusions."""

import argparse
from pathlib import Path
import shutil

import pandas as pd
import yaml

from genotype_samples import metadata_sample_ids


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('crosswalk', type=Path)
    args = parser.parse_args()
    root = Path('fastq/Amyg_Heart')
    config = yaml.safe_load(Path('config/config.yaml').read_text())
    crosswalk = pd.read_csv(args.crosswalk, sep='\t', dtype=str, keep_default_na=False)
    samples = pd.read_csv(root / 'sample_metadata.tsv', sep='\t', dtype=str, keep_default_na=False)
    joined = samples.merge(crosswalk, left_on='run_accession', right_on='ena_run_accession',
                           suffixes=('', '_crosswalk'), how='outer', validate='one_to_one', indicator=True)
    assert joined._merge.eq('both').all(), 'Crosswalk and ENA run sets differ'
    assert joined.archive_sample.eq(joined.public_sample_id).all(), 'Archive sample identities differ'
    assert joined.tissue.map({'Amyg': 'amygdala', 'Heart': 'heart'}).eq(joined.tissue_crosswalk).all()
    assert joined.study_accession.map({'PRJEB60407': 'E-MTAB-12693',
                                      'PRJEB60349': 'E-MTAB-12701'}).eq(joined.biostudies_accession).all()
    assert joined.plink_fid.eq('HSrats').all()
    assert joined.plink_iid.str.fullmatch(r'[0-9]+').all()
    joined['rfid'] = 'Outbred' + joined.plink_iid
    assert not joined.duplicated(['tissue', 'rfid']).any()
    settings = config['genotypes']
    candidates = metadata_sample_ids(settings['metadata'])
    first = None
    for chrom in settings['chromosomes']:
        psam = pd.read_csv(settings['source_prefix'].format(chrom=chrom) + '.psam',
                           sep=r'\s+', dtype=str).rename(columns={'#IID': 'IID'})
        if first is None:
            first = psam
        else:
            assert first.equals(psam), f'Chromosome {chrom} PSAM differs'
    available = first.merge(candidates, on='IID', how='left', validate='one_to_one')
    assert available.rfid.notna().all()
    counts = available.groupby('rfid').size()
    joined['genotype_count'] = joined.rfid.map(counts).fillna(0).astype(int)
    joined['genotype_status'] = joined.genotype_count.map({0: 'missing', 1: 'present'}).fillna('ambiguous')
    joined['in_genotype_metadata'] = joined.rfid.isin(candidates.rfid)
    joined = joined.drop(columns='_merge').sort_values(['tissue', 'rfid'])
    joined.to_csv(root / 'genotype_coverage.tsv', sep='\t', index=False)
    for tissue in ['Amyg', 'Heart']:
        rows = joined.loc[joined.tissue.eq(tissue)]
        ids = ''.join(r + '\n' for r in rows.rfid)
        maps = ''.join('\t'.join([r.fastq_1] + ([r.fastq_2] if r.fastq_2 else []) + [r.rfid]) + '\n'
                       for r in rows.itertuples())
        release = Path(config['version']) / tissue
        release.mkdir(parents=True, exist_ok=True)
        for directory, name, text in [
            (Path('samples/rat_ids_current'), f'rat_ids-{tissue}.txt', ids),
            (Path('samples/fastq_map_current'), f'fastq_map-{tissue}.txt', maps),
            (release, 'rat_ids.txt', ids), (release, 'fastq_map.txt', maps),
        ]:
            (directory / name).write_text(text)
        rows.to_csv(release / 'sample_metadata.tsv', sep='\t', index=False)
        print(tissue, rows.genotype_status.value_counts().to_dict())
    shutil.copyfile(args.crosswalk, root / 'metadata' / args.crosswalk.name)


if __name__ == '__main__':
    main()
