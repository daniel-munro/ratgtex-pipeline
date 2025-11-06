rule vcf_to_plink:
    """Get SNPs that are not monomorphic in a given set of samples."""
    input:
        vcf = lambda w: f"geno/{config['tissues'][w.tissue]['geno_dataset']}.vcf.gz",
        samples = "{version}/{tissue}/rat_ids.txt"
    output:
        multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam")
    params:
        prefix = "{version}/{tissue}/geno"
    shell:
        """
        plink2 --make-bed \
            --vcf {input.vcf} \
            --keep {input.samples} \
            --maf 0.01 \
            --mac 2 \
            --max-alleles 2 \
            --output-chr chrM \
            --out {params.prefix}
        """


rule tensorqtl_cis_nominal:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/phenos/output/expression.bed.gz",
        bedi = "{version}/{tissue}/phenos/output/expression.bed.gz.tbi",
        covar = "{version}/{tissue}/pheast/intermediate/covar/expression.covar.tsv"
    output:
        expand("{{version}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.chr{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}/nominal",
        out_prefix = "{tissue}",
    resources:
        mem_mb = 32000,
        runtime = '12h',
    shell:
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
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/phenos/output/expression.bed.gz",
        bedi = "{version}/{tissue}/phenos/output/expression.bed.gz.tbi",
        covar = "{version}/{tissue}/pheast/intermediate/covar/expression.covar.tsv"
    output:
        "{version}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}",
        out_prefix = "{tissue}",
    resources:
        mem_mb = lambda w, attempt: 32000 * 2**(attempt-1),
        runtime = '12h',
    retries: 2
    shell:
        # batch_size set due to "RuntimeError: CUDA out of memory"
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
        qtl = "{version}/{tissue}/pheast/output/qtl/expression.cis_qtl.txt.gz",
        nom = expand("{{version}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.chr{chrn}.parquet", chrn=range(1, 21))
    output:
        "{version}/{tissue}/{tissue}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "{version}/{tissue}/nominal/{tissue}"
    shell:
        "python3 scripts/tensorqtl_all_signif.py {input.qtl} {params.nom_prefix} {output} --fdr 0.05"


rule tensorqtl_all_cis_pvals:
    """Extract p-values for all tested cis-window SNPs per gene."""
    input:
        expand("{{version}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.chr{chrn}.parquet", chrn=range(1, 21))
    output:
        "{version}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz"
    params:
        nom_dir = "{version}/{tissue}/nominal"
    shell:
        "python3 scripts/tensorqtl_all_cis_pvals.py {params.nom_dir} {output}"


rule assemble_log_expression:
    """Pantry already assembles the TPM expression phenotype table, but this
    generates log2(readcount+1) tables for aFC and download.
    """
    input:
        kallisto = lambda w: expand("{{version}}/{{tissue}}/phenos/intermediate/expression/{rat_id}/abundance.tsv", rat_id=ids(w.tissue)),
        samples = "{version}/{tissue}/rat_ids.txt",
        ref_anno = f"{ANNO_PREFIX}.gtf",
    output:
        multiext("{version}/{tissue}/{tissue}.expr.log2.bed", ".gz", ".gz.tbi"),
    params:
        script_path = "{version}/{tissue}/phenos/scripts/assemble_bed.py",
        expr_dir = "{version}/{tissue}/phenos/intermediate/expression",
        bed = "{version}/{tissue}/{tissue}.expr.log2.bed",
        bed_iso = "{version}/{tissue}/{tissue}.iso.log2.bed",
    resources:
        mem_mb = 32000,
    shell:
        """
        python3 {params.script_path} expression \
            --input-dir {params.expr_dir} \
            --samples {input.samples} \
            --ref-anno {input.ref_anno} \
            --output-isoforms {params.bed_iso} \
            --output-expression {params.bed} \
            --units est_counts \
            --log2-expr
        rm {params.bed_iso}
        bgzip {params.bed}
        tabix {params.bed}.gz
        """


rule aFC:
    """Get effect size (allelic fold change) for top association per gene and all significant eQTLs."""
    input:
        vcf = lambda w: f"geno/{config['tissues'][w.tissue]['geno_dataset']}.vcf.gz",
        vcfi = lambda w: f"geno/{config['tissues'][w.tissue]['geno_dataset']}.vcf.gz.tbi",
        bed = "{version}/{tissue}/{tissue}.expr.log2.bed.gz",
        bedi = "{version}/{tissue}/{tissue}.expr.log2.bed.gz.tbi",
        qtl = "{version}/{tissue}/pheast/output/qtl/expression.cis_qtl.txt.gz",
        qtl_indep = "{version}/{tissue}/pheast/output/qtl/expression.cis_independent_qtl.txt.gz",
        covar = "{version}/{tissue}/pheast/intermediate/covar/expression.covar.tsv",
    output:
        "{version}/{tissue}/{tissue}.aFC.txt"
    resources:
        runtime = '12h'
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


