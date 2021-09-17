def ids(tissue):
    """Load list of rat IDs for a tissue dataset."""
    with open(f"{tissue}/rat_ids.txt", "r") as f:
        return f.read().splitlines()


# These steps are short and will not be submitted as cluster jobs:
localrules:
    individual_vcf,
    observed_snp_list,
    vcf_chr_list,
    index_vcf,
    index_bam,
    plink_sample_list,
    vcf_to_plink,
    expression_tsv_to_bed,
    expression_pcs,
    combine_covariates,
    empty_covariates,
    covariates,


rule all:
    """List any files here that you want to generate by default
    (i.e. without having to specify on command line when running snakemake).
    """
    input:
        # expand("Eye/markdup_out/{id}.bam", id=ids["Eye"][:10])
        # "Eye/Eye.cis_qtl_signif.txt.gz",
        # "Eye/Eye.cis_qtl_all_pvals.txt.gz",
        # "Eye/Eye.aFC.txt",
        # "Eye/nominal/Eye.cis_qtl_pairs.1.parquet",
        # "Eye/Eye.trans_qtl_pairs.txt.gz",


# TODO: add STAR step

# TODO: add RSEM step


rule sim_to_founders:
    """Calculate genetic similarity of each rat to each founder strain to use as covariates.
    Since all tissues share a VCF, this is calculated only once.
    """
    input:
        pop_vcf = "geno/ratgtex.vcf.gz",
        founder_vcf = "geno/founders.vcf.gz"
    output:
        "geno/sim_to_founders.txt"
    shell:
        "Rscript src/sim_to_founders.R {input.pop_vcf} {input.founder_vcf} {output}"


rule assemble_expression:
    """Combine RSEM output from all samples into log-count and TPM expression tables.
    Also computes inverse-quantile normalized values to use for eQTL mapping. The filtered
    version is to avoid a tensorQTL error on phenotypes with 1 nonzero value.
    """
    input:
        rsem = lambda w: expand("{{tissue}}/rsem_out/{rat_id}.genes.results.gz", rat_id=ids(w.tissue)),
        anno = "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        multiext("{tissue}/{tissue}.expr.log2.bed", ".gz", ".gz.tbi"),
        multiext("{tissue}/{tissue}.expr.tpm.bed", ".gz", ".gz.tbi"),
        multiext("{tissue}/{tissue}.expr.iqn.bed", ".gz", ".gz.tbi"),
        multiext("{tissue}/{tissue}.expr.iqn.filtered.bed", ".gz", ".gz.tbi")
    params:
        rsem_dir = "{tissue}/rsem_out",
        prefix = "{tissue}/{tissue}.expr"
    shell:
        """
        python3 src/assemble_expression.py {params.rsem_dir} {input.anno} {params.prefix}
        bgzip {params.prefix}.log2.bed
        bgzip {params.prefix}.tpm.bed
        bgzip {params.prefix}.iqn.bed
        bgzip {params.prefix}.iqn.filtered.bed
        tabix {params.prefix}.log2.bed.gz
        tabix {params.prefix}.tpm.bed.gz
        tabix {params.prefix}.iqn.bed.gz
        tabix {params.prefix}.iqn.filtered.bed.gz
        """


rule covariates:
    """Compute expression PCs and combine with genotype covariates."""
    input:
        geno = "geno/sim_to_founders.txt",
        bed = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
    output:
        "{tissue}/covar.txt"
    params:
        n_pcs = 20
    shell:
        "Rscript src/covariates.R {input.geno} {input.bed} {params.n_pcs} {output}"


rule vcf_to_plink:
    """Get SNPs that are not monomorphic in a given set of samples."""
    input:
        vcf = "geno/ratgtex.vcf.gz",
        samples = "{tissue}/rat_ids.txt"
    output:
        multiext("{tissue}/geno", ".bed", ".bim", ".fam")
    params:
        prefix = "{tissue}/geno"
    shell:
        """
        plink2 --make-bed \
        --vcf {input.vcf} \
        --keep {input.samples} \
        --maf 0.01 \
        --max-maf 0.99 \
        --out {params.prefix}
        """


rule tensorqtl_perm:
    """Map cis-eQTLs, determining significance using permutations.
    Outputs the top association per gene.
    """
    input:
        geno = multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{tissue}/covar.txt"
    output:
        "{tissue}/{tissue}.cis_qtl.txt.gz"
    params:
        geno_prefix = "{tissue}/geno",
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        """
        module load cuda
        # 
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --mode cis
        """


rule tensorqtl_independent:
    input:
        geno = multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{tissue}/covar.txt",
        cis = "{tissue}/{tissue}.cis_qtl.txt.gz"
    output:
        "{tissue}/{tissue}.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = "{tissue}/geno",
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        """
        module load cuda
        python3 src/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --mode cis_independent
        """


rule tensorqtl_nominal:
    input:
        geno = multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{tissue}/covar.txt"
    output:
        expand("{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = "{tissue}/geno",
        outdir = "{tissue}/nominal"
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        """
        module load cuda
        # pip install -e ~/tools/tensorqtl
        mkdir -p {params.outdir}
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.tissue} \
            --covariates {input.covar} \
            --output_dir {params.outdir} \
            --mode cis_nominal
        """


rule tensorqtl_trans:
    input:
        geno = multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{tissue}/covar.txt"
    output:
        "{tissue}/{tissue}.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "{tissue}/geno",
        outdir = "{tissue}"
    resources:
        walltime = 12,
        partition = "--partition=gpu",
    shell:
        # batch_size set due to "RuntimeError: CUDA out of memory"
        """
        module load cuda
        python3 -m tensorqtl \
            {params.geno_prefix} \
            {input.bed} \
            {wildcards.tissue} \
            --covariates {input.covar} \
            --output_dir {params.outdir} \
            --output_text \
            --batch_size 10000 \
            --mode trans
        """


rule tensorqtl_all_signif:
    input:
        perm = "{tissue}/{tissue}.cis_qtl.txt.gz",
        nom = expand("{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "{tissue}/{tissue}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "{tissue}/nominal/{tissue}"
    shell:
        "python3 src/tensorqtl_all_signif.py {input.perm} {params.nom_prefix} {output}"


rule tensorqtl_all_cis_pvals:
    input:
        expand("{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz"
    params:
        nom_dir = "{tissue}/nominal"
    shell:
        "python3 src/tensorqtl_all_cis_pvals.py {params.nom_dir} {output}"


rule aFC:
    input:
        vcf = "geno/ratgtex.vcf.gz",
        vcfi = "geno/ratgtex.vcf.gz.tbi",
        bed = "{tissue}/{tissue}.expr.log2.bed.gz",
        bedi = "{tissue}/{tissue}.expr.log2.bed.gz.tbi",
        qtl = "{tissue}/{tissue}.cis_qtl.txt.gz",
        qtl_indep = "{tissue}/{tissue}.cis_independent_qtl.txt.gz",
        covar = "{tissue}/covar.txt"
    output:
        "{tissue}/{tissue}.aFC.txt"
    resources:
        walltime = 12
    shell:
        """
        python3 ~/tools/aFC/aFC.py \
        --vcf {input.vcf} \
        --pheno {input.bed} \
        --qtl <(python3 src/prepare_qtl_for_afc.py {input.qtl} {input.qtl_indep}) \
        --cov {input.covar} \
        --log_xform 1 \
        --output {output}
        """

