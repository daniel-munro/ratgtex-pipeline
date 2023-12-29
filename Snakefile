from pathlib import Path


def ids(tissue):
    """Load list of rat IDs for a tissue dataset."""
    with open(f"{RN}/{tissue}/rat_ids.txt", "r") as f:
        return f.read().splitlines()


RN = "rn7"
if RN == "rn6":
    GENOME_PREFIX = "ref_rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel"
    ANNO_PREFIX = "ref_rn6/Rattus_norvegicus.Rnor_6.0.99"
elif RN == "rn7":
    GENOME_PREFIX = "ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel"
    ANNO_PREFIX = "ref_rn7/Rattus_norvegicus.mRatBN7.2.108"

configfile: 'config.yaml'
# For each tissue, specify: read_length, fastq_path, paired_end, geno_dataset

# These steps are short and will not be submitted as cluster jobs:
localrules:
    # index_vcf,
    prune_for_covar,
    covariates,

# Alignment, quality control, and expression steps can be found in these files:
include: "align.smk"
include: "qc.smk"
include: "expression.smk"
include: "splicing.smk"

# TISSUE = "IL"
# TISSUES = ["IL", "LHb", "NAcc", "OFC", "PL"]
# TISSUES = ["BLA", "NAcc2", "PL2"]
# TISSUES = ["Brain"]
# TISSUES = ["Adipose", "Liver"]
TISSUES = ["IL", "LHb", "NAcc", "OFC", "PL", "Eye", "BLA", "NAcc2", "PL2", "Adipose", "Liver", "Brain"]
rule all:
    """List any files here that you want to generate by default
    (i.e. without having to specify on command line when running snakemake).
    """
    input:
        # f"{RN}/{TISSUE}/qc/rna_to_geno_summary.tsv",
        # f"{RN}/{TISSUE}/qc/all_rats_summary.tsv",
        # f"{RN}/{TISSUE}/{TISSUE}.cis_qtl_signif.txt.gz",
        # f"{RN}/{TISSUE}/{TISSUE}.cis_qtl_all_pvals.txt.gz",
        # f"{RN}/{TISSUE}/{TISSUE}.aFC.txt",
        # f"{RN}/{TISSUE}/{TISSUE}.trans_qtl_pairs.txt.gz",
        # f"{RN}/{TISSUE}/splice/{TISSUE}_splice.cis_independent_qtl.txt.gz",
        # expand("{rn}/{tissue}/qc/rna_to_geno_summary.tsv", rn=RN, tissue=TISSUES),
        # expand("{rn}/{tissue}/qc/all_rats_summary.tsv", rn=RN, tissue=TISSUES),
        # expand("{rn}/{tissue}/qc/{tissue}.sex_concordance.txt", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.expr.tpm.bed.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.cis_qtl.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.cis_independent_qtl.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.cis_qtl_signif.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.aFC.txt", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz", rn=RN, tissue=TISSUES),
        expand("{rn}/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz", rn=RN, tissue=TISSUES),


rule vcf_to_plink:
    """Get SNPs that are not monomorphic in a given set of samples."""
    input:
        vcf = lambda w: f"geno_{RN}/{config[w.tissue]['geno_dataset']}.vcf.gz",
        samples = "{rn}/{tissue}/rat_ids.txt"
    output:
        multiext("{rn}/{tissue}/geno", ".bed", ".bim", ".fam")
    params:
        prefix = "{rn}/{tissue}/geno"
    shell:
        """
        plink2 --make-bed \
            --vcf {input.vcf} \
            --keep {input.samples} \
            --maf 0.01 \
            --mac 2 \
            --max-alleles 2 \
            --out {params.prefix}
        """


rule prune_for_covar:
    """Prune genotypes to compute covariate PCs.
    --indep-pairwise parameters are based on GTEx methods.
    """
    input:
        multiext("{rn}/{tissue}/geno", ".bed", ".bim", ".fam")
    output:
        "{rn}/{tissue}/covar/geno.vcf.gz"
    params:
        prefix = "{rn}/{tissue}/geno",
        pruned_dir = "{rn}/{tissue}/covar",
        pruned_prefix = "{rn}/{tissue}/covar/geno"
    shell:
        # --geno 0.05 filters variants with >5% missing values (the rest will be imputed)
        # Updated to --geno 0.00 due to cases where one or a few samples have many fewer genotyped
        # SNPs (e.g. using different genotyping methods). Without subsetting to shared SNPs, those
        # samples couldn't be imputed later due to having so many missing values.
        """
        mkdir -p {params.pruned_dir}
        plink2 \
            --bfile {params.prefix} \
            --geno 0.00 \
            --maf 0.05 \
            --indep-pairwise 200 100 0.1 \
            --out {params.pruned_prefix}
        plink2 \
            --bfile {params.prefix} \
            --extract {params.pruned_prefix}.prune.in \
            --export vcf bgz id-paste=iid \
            --out {params.pruned_prefix}
        """


rule covariates:
    """Compute genotype and expression PCs and combine."""
    input:
        vcf = "{rn}/{tissue}/covar/geno.vcf.gz",
        bed = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
    output:
        "{rn}/{tissue}/covar.txt"
    params:
        n_geno_pcs = 5,
        n_expr_pcs = 20
    shell:
        "Rscript scripts/covariates.R {input.vcf} {input.bed} {params.n_geno_pcs} {params.n_expr_pcs} {output}"


rule tensorqtl_perm:
    """Map cis-eQTLs, determining significance using permutations.
    Outputs the top association per gene.
    """
    input:
        geno = multiext("{rn}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{rn}/{tissue}/covar.txt"
    output:
        "{rn}/{tissue}/{tissue}.cis_qtl.txt.gz"
    params:
        geno_prefix = "{rn}/{tissue}/geno",
    resources:
        walltime = 20,
        partition = "--partition=gpu",
    # retries: 5
    shell:
        """
        module load cuda
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --mode cis
        """


rule tensorqtl_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-eQTLs per gene."""
    input:
        geno = multiext("{rn}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{rn}/{tissue}/covar.txt",
        cis = "{rn}/{tissue}/{tissue}.cis_qtl.txt.gz"
    output:
        "{rn}/{tissue}/{tissue}.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = "{rn}/{tissue}/geno",
    resources:
        walltime = 20,
        partition = "--partition=gpu",
    # retries: 5
    shell:
        """
        module load cuda
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --mode cis_independent
        """


rule tensorqtl_nominal:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        geno = multiext("{rn}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{rn}/{tissue}/covar.txt"
    output:
        expand("{{rn}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = "{rn}/{tissue}/geno",
        outdir = "{rn}/{tissue}/nominal",
        out_prefix = "{tissue}",
    resources:
        walltime = 12,
        # partition = "--partition=gpu",
    shell:
        # module load cuda
        """
        mkdir -p {params.outdir}
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            {params.out_prefix} \
            --covariates {input.covar} \
            --output_dir {params.outdir} \
            --mode cis_nominal
        """


rule tensorqtl_trans:
    """Map trans-eQTLs (without significance testing)."""
    input:
        geno = multiext("{rn}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{rn}/{tissue}/covar.txt"
    output:
        "{rn}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "{rn}/{tissue}/geno",
        outdir = "{rn}/{tissue}",
        out_prefix = "{tissue}",
    resources:
        walltime = 12,
        # partition = "--partition=gpu",
    shell:
        # batch_size set due to "RuntimeError: CUDA out of memory"
        # module load cuda
        """
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            {params.out_prefix} \
            --covariates {input.covar} \
            --output_dir {params.outdir} \
            --output_text \
            --batch_size 10000 \
            --mode trans
        """


rule tensorqtl_all_signif:
    """Extract all significant cis SNP-gene pairs."""
    input:
        perm = "{rn}/{tissue}/{tissue}.cis_qtl.txt.gz",
        nom = expand("{{rn}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "{rn}/{tissue}/{tissue}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "{rn}/{tissue}/nominal/{tissue}"
    shell:
        "python3 scripts/tensorqtl_all_signif.py {input.perm} {params.nom_prefix} {output} --fdr 0.05"


rule tensorqtl_all_cis_pvals:
    """Extract p-values for all tested cis-window SNPs per gene."""
    input:
        expand("{{rn}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "{rn}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz"
    params:
        nom_dir = "{rn}/{tissue}/nominal"
    shell:
        "python3 scripts/tensorqtl_all_cis_pvals.py {params.nom_dir} {output}"


rule aFC:
    """Get effect size (allelic fold change) for top association per gene and all significant eQTLs."""
    input:
        vcf = lambda w: f"geno_{RN}/{config[w.tissue]['geno_dataset']}.vcf.gz",
        vcfi = lambda w: f"geno_{RN}/{config[w.tissue]['geno_dataset']}.vcf.gz.tbi",
        bed = "{rn}/{tissue}/{tissue}.expr.log2.bed.gz",
        bedi = "{rn}/{tissue}/{tissue}.expr.log2.bed.gz.tbi",
        qtl = "{rn}/{tissue}/{tissue}.cis_qtl.txt.gz",
        qtl_indep = "{rn}/{tissue}/{tissue}.cis_independent_qtl.txt.gz",
        covar = "{rn}/{tissue}/covar.txt"
    output:
        "{rn}/{tissue}/{tissue}.aFC.txt"
    resources:
        walltime = 12
    shell:
        """
        python3 tools/aFC/aFC.py \
            --vcf {input.vcf} \
            --pheno {input.bed} \
            --qtl <(python3 scripts/prepare_qtl_for_afc.py {input.qtl} {input.qtl_indep}) \
            --cov {input.covar} \
            --log_xform 1 \
            --output {output}
        """


# rule index_vcf:
#     """Generate an index file for a VCF file."""
#     input:
#         "{base}.vcf.gz"
#     output:
#         "{base}.vcf.gz.tbi"
#     shell:
#         "tabix -p vcf {input}"


