"""Plot STAR sample quality for configured v5 datasets without changing membership."""

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

project_root = (Path(snakemake.scriptdir).resolve().parents[1]
                if 'snakemake' in globals() else Path(__file__).resolve().parents[2])
sys.path.insert(0, str(project_root / 'scripts' / 'qc'))
from stats_from_star_logs import parse_star_log

if 'snakemake' in globals():
    config = snakemake.config
    out_dir = Path(snakemake.output.metrics).parent
else:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=Path('config/config.yaml'))
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    with args.config.open() as f:
        config = yaml.safe_load(f)
    out_dir = args.out_dir
records = []
for tissue, datasets in config['tissues'].items():
    for dataset in datasets:
        base = Path(config['version']) / dataset
        for sample in (base / 'rat_ids.txt').read_text().splitlines():
            stats = parse_star_log(base / 'star_out' / f'{sample}.Log.final.out')
            records.append(dict(
                tissue=tissue, dataset=dataset, sample_id=sample,
                input_reads=int(stats['Number of input reads']),
                uniquely_mapped_reads=int(stats['Uniquely mapped reads number']),
                uniquely_mapped_percent=float(stats['Uniquely mapped reads %'].rstrip('%')),
            ))
df = pd.DataFrame(records)
out_dir.mkdir(parents=True, exist_ok=True)
df.to_csv(out_dir / 'sample_metrics.tsv', sep='\t', index=False)
thresholds = [500_000, 1_000_000, 2_000_000]
summary = df.groupby('tissue', sort=False).size().rename('samples').to_frame()
for threshold in thresholds:
    column = f'below_{threshold}'
    df[column] = df.uniquely_mapped_reads < threshold
    summary[column] = df.groupby('tissue', sort=False)[column].sum()
summary.to_csv(out_dir / 'candidate_cutoff_counts.tsv', sep='\t')
df.loc[df.below_2000000].sort_values(['tissue', 'uniquely_mapped_reads']).to_csv(
    out_dir / 'samples_below_2million.tsv', sep='\t', index=False)

tissues = list(config['tissues'])
ncols = 4
nrows = int(np.ceil(len(tissues) / ncols))
for metric, title, xlabel, filename in [
    ('uniquely_mapped_reads', 'Uniquely mapped depth by tissue',
     'Uniquely mapped reads / read pairs (log scale)', 'uniquely_mapped_depth'),
    ('uniquely_mapped_percent', 'Unique mapping percentage by tissue',
     'Uniquely mapped (%)', 'unique_mapping_percent'),
]:
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 3 * nrows), squeeze=False)
    if metric == 'uniquely_mapped_reads':
        if (df[metric] <= 0).any():
            raise ValueError('Log-depth plot requires positive uniquely mapped counts')
        bins = np.geomspace(min(df[metric].min(), 100_000) * .9,
                            df[metric].max() * 1.1, 36)
    else:
        bins = np.linspace(0, 100, 31)
    for ax, tissue in zip(axes.flat, tissues):
        sub = df.loc[df.tissue == tissue]
        groups = list(sub.groupby('dataset', sort=False))
        ax.hist([g[metric].to_numpy() for _, g in groups], bins=bins,
                stacked=True, label=[d for d, _ in groups], edgecolor='white', linewidth=.3)
        if metric == 'uniquely_mapped_reads':
            ax.set_xscale('log')
            for cutoff, color in zip(thresholds, ['#e69f00', '#d55e00', '#777777']):
                ax.axvline(cutoff, color=color, linestyle='--', linewidth=1,
                           label=f'{cutoff/1e6:g}M')
            ax.set_xlim(bins[0], bins[-1])
            counts = summary.loc[tissue]
            ax.text(.02, .95, '\n'.join(f'<{t/1e6:g}M: {int(counts[f"below_{t}"])}'
                                      for t in thresholds),
                    transform=ax.transAxes, va='top', fontsize=8)
        else:
            ax.set_xlim(0, 100)
        ax.set_title(f'{tissue} (n={len(sub)})')
        ax.set_ylabel('Samples')
        ax.set_xlabel(xlabel, fontsize=8)
        ax.legend(fontsize=7, loc='upper right')
    for ax in list(axes.flat)[len(tissues):]:
        ax.set_visible(False)
    fig.suptitle(title + '\nCurrent v5 membership; candidate cutoffs only, no exclusions applied', fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, .95))
    for ext in ['png', 'pdf']:
        fig.savefig(out_dir / f'{filename}.{ext}', dpi=170)
    plt.close(fig)
print(summary.to_string())
print(f'\nTotal samples: {len(df)}')
print(df.loc[df.below_1000000, ['tissue', 'dataset', 'sample_id', 'uniquely_mapped_reads',
                             'uniquely_mapped_percent']].to_string(index=False))
