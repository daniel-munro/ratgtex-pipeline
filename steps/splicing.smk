localrules:
    exon_table,


rule regtools_junctions:
    """Much faster than LeafCutter's version, so even they now recommend using regtools."""
    input:
        bam = "{version}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bai = "{version}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "{version}/{tissue}/splice/junc/{rat_id}.junc.gz"
    params:
        junc_dir = "{version}/{tissue}/splice/junc",
        intermediate = "{version}/{tissue}/splice/junc/{rat_id}.junc"
    shell:
        """
        mkdir -p {params.junc_dir}
        regtools junctions extract \
            -a 8 \
            -m 50 \
            -M 500000 \
            -s XS \
            -o {params.intermediate} \
            {input.bam}
        gzip {params.intermediate}
        """


rule exon_table:
    input:
        f"{ANNO_PREFIX}.genes.gtf"
    output:
        f"{ANNO_PREFIX}.genes.exons.tsv"
    shell:
        "Rscript scripts/splice/exon_table.R {input} {output}"


rule splice_bed:
    """Note: fails if intermediate files are already present.
    
    For some reason the phenotype groups file produced by
    cluster_prepare_fastqtl.py includes a line '\t0', so that is removed so
    tensorQTL nominal mode will work.
    """
    input:
        junc = lambda w: expand("{{version}}/{{tissue}}/splice/junc/{rat_id}.junc.gz", rat_id=ids(w.tissue)),
        exons = F"{ANNO_PREFIX}.genes.exons.tsv",
        gtf = F"{ANNO_PREFIX}.genes.gtf",
    output:
        bed = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        clust = "{version}/{tissue}/splice/{tissue}.leafcutter.clusters_to_genes.txt",
        pcs = "{version}/{tissue}/splice/{tissue}.leafcutter.PCs.txt",
        groups = "{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
    params:
        prefix = "{tissue}",
        script_dir = "scripts/splice",
        tmpdir = "{version}/{tissue}/splice/clust",
    shell:
        """
        mkdir -p {params.tmpdir}
        echo {input.junc} | tr ' ' '\n' > {params.tmpdir}/juncfiles.txt
        python3 scripts/splice/cluster_prepare_fastqtl.py \
            {params.tmpdir}/juncfiles.txt \
            {input.exons} \
            {input.gtf} \
            {params.prefix} \
            --leafcutter_dir {params.script_dir} \
            --output_dir {params.tmpdir}
        mv {params.tmpdir}/{wildcards.tissue}.leafcutter.phenotype_groups.txt {params.tmpdir}/{wildcards.tissue}.leafcutter.phenotype_groups.tmp.txt
        grep -P -v '^\t0' {params.tmpdir}/{wildcards.tissue}.leafcutter.phenotype_groups.tmp.txt > {params.tmpdir}/{wildcards.tissue}.leafcutter.phenotype_groups.txt
        rm {params.tmpdir}/{wildcards.tissue}.leafcutter.phenotype_groups.tmp.txt
        mv -i {params.tmpdir}/{wildcards.tissue}.leafcutter* {params.tmpdir}/..
        rm -r {params.tmpdir}
        """


rule splice_covariates:
    """Compute genotype and splicing PCs and combine."""
    input:
        geno = multiext("{version}/{tissue}/covar/geno_pruned", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz",
    output:
        covar = "{version}/{tissue}/splice/{tissue}.covar_splice.txt"
    params:
        pruned_prefix = "{version}/{tissue}/covar/geno_pruned",
        n_geno_pcs = 5,
        n_pheno_pcs = 10
    shell:
        "Rscript scripts/covariates.R {params.pruned_prefix} {input.bed} {params.n_geno_pcs} {params.n_pheno_pcs} {output.covar}"


rule tensorqtl_cis_splice:
    input:
        geno = lambda w: multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{version}/{tissue}/splice/{tissue}.covar_splice.txt",
        groups = "{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
    output:
        "{version}/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}/splice"
    resources:
        runtime = '20h',
    shell:
        """
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --groups {input.groups} \
            --mode cis
        """


rule tensorqtl_cis_independent_splice:
    input:
        geno = lambda w: multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{version}/{tissue}/splice/{tissue}.covar_splice.txt",
        groups = "{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
        cis = "{version}/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz",
    output:
        "{version}/{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}/splice"
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
            --groups {input.groups} \
            --mode cis_independent
        """


rule tensorqtl_trans_splice:
    """Map trans-sQTLs (without significance testing)."""
    input:
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{version}/{tissue}/splice/{tissue}.covar_splice.txt",
    output:
        "{version}/{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}/splice",
        out_prefix = "{tissue}_splice",
    resources:
        runtime = '12h',
    shell:
        # batch_size set due to "RuntimeError: CUDA out of memory"
        # tensorQTL doesn't use group info for trans mapping
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


rule tensorqtl_cis_nominal_splice:
    """Get summary statistics for all tested cis-window SNPs per gene."""
    input:
        geno = multiext("{version}/{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{version}/{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{version}/{tissue}/splice/{tissue}.covar_splice.txt",
        groups = "{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
    output:
        expand("{{version}}/{{tissue}}/splice/nominal/{{tissue}}_splice.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
    params:
        geno_prefix = "{version}/{tissue}/geno",
        outdir = "{version}/{tissue}/splice/nominal",
        out_prefix = "{tissue}_splice",
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
            --phenotype_groups {input.groups} \
            --output_dir {params.outdir} \
            --mode cis_nominal
        """


rule tensorqtl_all_signif_splice:
    """Extract all significant cis SNP-gene pairs."""
    input:
        perm = "{version}/{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz",
        nom = expand("{{version}}/{{tissue}}/splice/nominal/{{tissue}}_splice.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21)),
        groups = "{version}/{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
    output:
        "{version}/{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz",
    params:
        nom_prefix = "{version}/{tissue}/splice/nominal/{tissue}_splice",
    shell:
        """
        python3 scripts/tensorqtl_all_signif.py \
            {input.perm} \
            {params.nom_prefix} \
            {output} \
            --groups {input.groups} \
            --fdr 0.05
        """

