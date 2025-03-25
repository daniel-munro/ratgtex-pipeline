from pathlib import Path


def ids(tissue):
    """Load list of rat IDs for a tissue dataset."""
    with open(f"{VERSION}/{tissue}/rat_ids.txt", "r") as f:
        return f.read().splitlines()


configfile: 'config.yaml'
# For each tissue, specify: read_length, fastq_path, paired_end, geno_dataset

# These steps are short and will not be submitted as cluster jobs:
# (Each included snakefile specifies its own localrules)
localrules:

# Alignment, quality control, and expression steps can be found in these files:
include: "steps/setup.smk"
include: "steps/align.smk"
include: "steps/qc.smk"
include: "steps/expression.smk"
include: "steps/splicing.smk"
include: "steps/qtl.smk"

rule all:
    """List any files here that you want to generate by default
    (i.e. without having to specify on command line when running snakemake).
    """
    input:
        expand("{v}/{tissue}/qc/rna_to_geno_summary.tsv", v=VERSION, tissue=TISSUES),
        # expand("{v}/{tissue}/qc/all_rats_summary.tsv", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/qc/{tissue}.sex_concordance.txt", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.expr.tpm.bed.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.cis_qtl.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.cis_independent_qtl.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.cis_qtl_signif.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.aFC.txt", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz", v=VERSION, tissue=TISSUES),
        expand("{v}/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz", v=VERSION, tissue=TISSUES),


