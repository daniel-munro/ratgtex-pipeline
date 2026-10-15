"""Shared interpretation of RatGTEx stage configuration."""

def tissues(config):
    return config['tissues']


def datasets(config):
    return list(dict.fromkeys(d for members in tissues(config).values() for d in members))


def write_config_lists(config, directory):
    """Export selected dataset and tissue IDs beside the main configuration."""
    for filename, names in [('datasets.txt', datasets(config)), ('tissues.txt', tissues(config))]:
        (directory / filename).write_text(''.join(f'{name}\n' for name in sorted(names, key=str.casefold)))


def phenotype_config(config, name, root):
    members = tissues(config).get(name, [name])
    settings = config['datasets'][members[0]]
    project = root / config['version'] / name
    result = dict(config['pantry']['phenotyping'])
    result.update(read_length=settings['read_length'], fastq_dir='/',
                  fastq_map=str(project / 'fastq_map.txt'),
                  samples_file=str(project / 'rat_ids.txt'),
                  ref_genome=str(root / config['ref_genome']),
                  ref_anno=str(root / config['ref_anno']),
                  intermediate_ref_dir=str(root / 'ref'))
    # Pairing is inferred per FASTQ-map row, including mixed-layout tissues.
    producer_options = result.pop('modality_groups', {})
    result['modality_groups'] = {}
    for modality, spec in config['modalities'].items():
        files = [f'{modality}.bed.gz', f'{modality}.bed.gz.tbi']
        if spec['grouped']:
            files.append(f'{modality}.phenotype_groups.txt')
        producer = spec['producer']
        result['modality_groups'].setdefault(
            producer, dict(producer_options.get(producer, {}), files=[]))['files'].extend(files)
    return result


def pheast_config(config, name, root, analysis):
    project = root / config['version'] / name
    result = dict(config['pantry']['pheast'])
    if 'twas_snps' in result:
        result['twas_snps'] = str(root / result['twas_snps'])
    result.update(phenotype_dir=str(project / 'phenos/output'),
                  geno_prefix=str(project / 'geno'), samples_file=str(project / 'rat_ids.txt'))
    result['modalities'] = {m: {'grouped': s['grouped']} for m, s in config['modalities'].items()}
    if analysis == 'qtl' and config['cross_modality']:
        result['modalities']['cross_modality'] = {'grouped': True}
    files = (['{modality}.cis_qtl.txt.gz', '{modality}.cis_independent_qtl.txt.gz']
             if analysis == 'qtl' else ['{modality}.tar.bz2'])
    result['analyses'] = {analysis: {'files': files}}
    return result
