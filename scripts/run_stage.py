#!/usr/bin/env python3
"""Set up Pantry projects and run RatGTEx steps before and after manual Pantry runs."""
import argparse
import json
import shlex
import shutil
import subprocess
from pathlib import Path

import yaml

from workflow_config import datasets, phenotype_config, pheast_config, write_config_lists

ROOT = Path(__file__).resolve().parents[1]
STAGES = ['setup', 'align', 'identity_qc', 'genotypes', 'cross_modality', 'qtl_extra']


def read_ids(path):
    ids = path.read_text().splitlines()
    if len(ids) != len(set(ids)):
        raise ValueError(f'Duplicate animal IDs in {path}')
    return ids


def link(source, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.symlink_to(source)


def setup(config, root):
    """Create new project copies only; never replace an existing Pantry project."""
    version = root / config['version']
    source = Path(config['pantry']['source']).expanduser().resolve()
    source_revision = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
    source_status = subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True)
    # Preflight destination collisions before making changes. Inputs are loaded normally.
    names = list(dict.fromkeys(datasets(config) + list(config['tissues'])))
    for name in names:
        for module in ['phenos', 'pheast']:
            if (version / name / module).exists():
                raise FileExistsError(f'{version / name / module}: setup creates new projects only')
    sample_ids = {d: read_ids(version / d / 'rat_ids.txt') for d in datasets(config)}
    maps = {}
    for d in datasets(config):
        fastq_dir = (root / config['datasets'][d]['fastq_path']).resolve()
        rows = []
        for line in (version / d / 'fastq_map.txt').read_text().splitlines():
            fields = line.split('\t')
            rows.append([str(fastq_dir / p) for p in fields[:-1]] + [fields[-1]])
        maps[d] = rows
    for tissue, members in config['tissues'].items():
        combined = [animal for d in members for animal in sample_ids[d]]
        if len(combined) != len(set(combined)):
            raise ValueError(f'{tissue}: overlapping animals across datasets require an explicit library-merging decision')
        if members != [tissue]:
            if tissue in config['datasets']:
                raise ValueError(f'{tissue} is both a source dataset and a distinct merged tissue')
            (version / tissue).mkdir(parents=True, exist_ok=True)
            (version / tissue / 'rat_ids.txt').write_text(''.join(f'{s}\n' for s in sorted(combined)))
            (version / tissue / 'fastq_map.txt').write_text(''.join(
                '\t'.join(row) + '\n' for d in members for row in maps[d]))
    for name in names:
        project = version / name
        for module, source_module in [('phenos', 'phenotyping'), ('pheast', 'pheast')]:
            shutil.copytree(source / source_module, project / module,
                            ignore=shutil.ignore_patterns('input', 'intermediate', 'output', '.snakemake', '__pycache__'))
        write_project_configs(config, name, root)
        shutil.copyfile(project / 'pheast/config.qtl.yml', project / 'pheast/config.yml')
        (project / 'ratgtex_setup.json').write_text(json.dumps(
            {'pantry_source': str(source), 'pantry_revision': source_revision,
             'pantry_status': source_status, 'ratgtex_config': config}, indent=2) + '\n')
        if name in sample_ids:
            for animal in sample_ids[name]:
                bam = project / f'star_out/{animal}.Aligned.sortedByCoord.out.bam'
                for suffix in ['', '.bai']:
                    link(Path(str(bam) + suffix), project / f'phenos/intermediate/bam/{animal}.bam{suffix}')
    for tissue, members in config['tissues'].items():
        if members == [tissue]:
            continue
        for d in members:
            for animal in sample_ids[d]:
                for producer in dict.fromkeys(s['producer'] for s in config['modalities'].values()):
                    # This mapping is an explicit interface to Pantry sample-level products.
                    for template in config['pantry']['merge_products'][producer]:
                        relative = template.format(sample=animal)
                        link(version / d / 'phenos/intermediate' / relative,
                             version / tissue / 'phenos/intermediate' / relative)


def write_project_configs(config, name, root):
    project = root / config['version'] / name

    def write_config(path, settings):
        path.write_text(yaml.safe_dump(settings, sort_keys=False))

    settings = phenotype_config(config, name, root)
    if name in datasets(config):
        settings['fastq_dir'] = str((root / config['datasets'][name]['fastq_path']).resolve())
    write_config(project / 'phenos/config.yml', settings)
    original = (project / 'pheast/Snakefile').read_text()
    if original.count("configfile: 'config.yml'") != 1:
        raise ValueError('Expected one configfile declaration in the copied Pantry Pheast Snakefile')
    for analysis in ['qtl', 'twas']:
        write_config(project / f'pheast/config.{analysis}.yml', pheast_config(config, name, root, analysis))
        (project / f'pheast/Snakefile.{analysis}').write_text(
            original.replace("configfile: 'config.yml'", f"configfile: 'config.{analysis}.yml'"))


def run(config, stage, config_path, args):
    if stage not in STAGES:
        raise ValueError(f'Unknown RatGTEx stage: {stage}. Run Pantry phenotyping and Pheast directly in their project directories.')
    if stage == 'setup':
        if args.dry_run or args.plan:
            print('Setup preview: export dataset/tissue lists and create new Pantry projects for ' + ', '.join(
                dict.fromkeys(datasets(config) + list(config['tissues']))))
        else:
            setup(config, ROOT)
            write_config_lists(config, config_path.parent)
    else:
        command = ['snakemake', stage, '--snakefile', 'Snakefile', '--profile', args.profile,
                   '--jobs', str(args.jobs), '--configfile', str(config_path)]
        if args.dry_run:
            command.append('--dry-run')
        print(f'[{ROOT}] {shlex.join(command)}', flush=True)
        if not args.plan:
            subprocess.run(command, cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=STAGES)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/config.yaml')
    parser.add_argument('--profile', default='slurm')
    parser.add_argument('--jobs', type=int, default=10)
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('--plan', action='store_true', help='Print commands without invoking Snakemake or writing files')
    args = parser.parse_args()
    path = args.config.resolve()
    config = yaml.safe_load(path.read_text())
    run(config, args.stage, path, args)


if __name__ == '__main__':
    main()
