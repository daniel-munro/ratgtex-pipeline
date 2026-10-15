rule combine_modalities:
    input:
        lambda w: expand('{v}/{tissue}/phenos/output/{modality}.bed.gz',
                         v=VERSION, tissue=w.tissue, modality=MODALITIES)
    output:
        bed = '{version}/{tissue}/phenos/output/cross_modality.bed.gz',
        index = '{version}/{tissue}/phenos/output/cross_modality.bed.gz.tbi',
        groups = '{version}/{tissue}/phenos/output/cross_modality.phenotype_groups.txt'
    params:
        project = '{version}/{tissue}',
        modalities = MODALITIES
    shell:
        'python3 scripts/combine_modalities.py {params.project:q} {params.modalities:q}'

rule cross_modality:
    input:
        expand('{v}/{tissue}/phenos/output/cross_modality.{suffix}',
               v=VERSION, tissue=TISSUES,
               suffix=['bed.gz', 'bed.gz.tbi', 'phenotype_groups.txt']) if config['cross_modality'] else []
