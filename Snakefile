from pathlib import Path


def ids(tissue):
    """Load list of rat IDs for a tissue dataset."""
    with open(f"{VERSION}/{tissue}/rat_ids.txt", "r") as f:
        return f.read().splitlines()


if not workflow.overwrite_configfiles:
    configfile: 'config/config.yaml'
# These steps are short and will not be submitted as cluster jobs:
# (Each included snakefile specifies its own localrules)
localrules:

# Alignment, quality control, and expression steps can be found in these files:
include: "steps/setup.smk"
include: "steps/genotypes.smk"
include: "steps/align.smk"
include: "steps/qc.smk"
include: "steps/qtl.smk"
include: "steps/phenotypes.smk"

def stage_inputs(stage):
    if stage == 'genotypes':
        return list(rules.genotypes.input)
    if stage == 'align':
        return expand('{v}/{tissue}/qc/{tissue}.star_stats.tsv', v=VERSION, tissue=DATASETS) + list(rules.qc_sample_quality_report.output)
    if stage == 'identity_qc':
        return expand('{v}/{tissue}/qc/rna_to_geno_summary.tsv', v=VERSION, tissue=DATASETS) + list(rules.qc_identity_report.output)
    if stage == 'qtl_extra':
        return (
            expand('{v}/{tissue}/{tissue}.{modality}.cis_qtl_signif.txt.gz', v=VERSION, tissue=TISSUES, modality=MODALITIES)
            + expand('{v}/{tissue}/{tissue}.expression.cis_qtl_all_pvals.tsv.gz', v=VERSION, tissue=TISSUES)
            + expand('{v}/{tissue}/{tissue}.aFC.tsv', v=VERSION, tissue=TISSUES)
            + expand('{v}/{tissue}/{tissue}.{modality}.trans_qtl_pairs.txt.gz', v=VERSION, tissue=TISSUES, modality=MODALITIES))
    raise ValueError(f'Unknown top-level stage: {stage}')


rule all:
    default_target: True
    input: stage_inputs(config['default_stage'])

rule align:
    input: stage_inputs('align')

rule identity_qc:
    input: stage_inputs('identity_qc')

rule qtl_extra:
    input: stage_inputs('qtl_extra')
