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
            --out {params.prefix}
        """


rule prune_for_covar:
    """Prune genotypes to compute covariate PCs.
    --indep-pairwise parameters are based on GTEx methods.
    """
    input:
        multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam")
    output:
        "{version}/{tissue}/covar/geno.vcf.gz"
    params:
        prefix = "{version}/{tissue}/geno",
        pruned_dir = "{version}/{tissue}/covar",
        pruned_prefix = "{version}/{tissue}/covar/geno"
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
        vcf = "{version}/{tissue}/covar/geno.vcf.gz",
        bed = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
    output:
        "{version}/{tissue}/covar.txt"
    params:
        n_geno_pcs = 5,
        n_expr_pcs = 20
    shell:
        "Rscript scripts/covariates.R {input.vcf} {input.bed} {params.n_geno_pcs} {params.n_expr_pcs} {output}"


rule tensorqtl_cis:
    """Map cis-eQTLs, determining significance using permutations.
    Outputs the top association per gene.
    """
    input:
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{version}/{tissue}/covar.txt"
    output:
        "{version}/{tissue}/{tissue}.cis_qtl.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
    resources:
        runtime = '12h',
    shell:
        """
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --mode cis
        """


rule tensorqtl_cis_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-eQTLs per gene."""
    input:
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{version}/{tissue}/covar.txt",
        cis = "{version}/{tissue}/{tissue}.cis_qtl.txt.gz"
    output:
        "{version}/{tissue}/{tissue}.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
    resources:
        runtime = '20h',
    shell:
        """
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            --mode cis_independent
        """


rule tensorqtl_cis_nominal:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{version}/{tissue}/covar.txt"
    output:
        expand("{{version}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}/nominal",
        out_prefix = "{tissue}",
    resources:
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
        bed = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz",
        bedi = "{version}/{tissue}/{tissue}.expr.iqn.filtered.bed.gz.tbi",
        covar = "{version}/{tissue}/covar.txt"
    output:
        "{version}/{tissue}/{tissue}.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}",
        out_prefix = "{tissue}",
    resources:
        runtime = '12h',
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
        perm = "{version}/{tissue}/{tissue}.cis_qtl.txt.gz",
        nom = expand("{{version}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "{version}/{tissue}/{tissue}.cis_qtl_signif.txt.gz"
    params:
        nom_prefix = "{version}/{tissue}/nominal/{tissue}"
    shell:
        "python3 scripts/tensorqtl_all_signif.py {input.perm} {params.nom_prefix} {output} --fdr 0.05"


rule tensorqtl_all_cis_pvals:
    """Extract p-values for all tested cis-window SNPs per gene."""
    input:
        expand("{{version}}/{{tissue}}/nominal/{{tissue}}.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    output:
        "{version}/{tissue}/{tissue}.cis_qtl_all_pvals.txt.gz"
    params:
        nom_dir = "{version}/{tissue}/nominal"
    shell:
        "python3 scripts/tensorqtl_all_cis_pvals.py {params.nom_dir} {output}"


rule aFC:
    """Get effect size (allelic fold change) for top association per gene and all significant eQTLs."""
    input:
        vcf = lambda w: f"geno/{config['tissues'][w.tissue]['geno_dataset']}.vcf.gz",
        vcfi = lambda w: f"geno/{config['tissues'][w.tissue]['geno_dataset']}.vcf.gz.tbi",
        bed = "{version}/{tissue}/{tissue}.expr.log2.bed.gz",
        bedi = "{version}/{tissue}/{tissue}.expr.log2.bed.gz.tbi",
        qtl = "{version}/{tissue}/{tissue}.cis_qtl.txt.gz",
        qtl_indep = "{version}/{tissue}/{tissue}.cis_independent_qtl.txt.gz",
        covar = "{version}/{tissue}/covar.txt"
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


