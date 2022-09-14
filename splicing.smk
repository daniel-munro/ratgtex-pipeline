localrules:
    exon_table,
    splice_covariates,


rule regtools_junctions:
    """Much faster than LeafCutter's version, so even they now recommend using regtools."""
    input:
        bam = "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bai = "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam.bai"
    output:
        "{tissue}/splice/junc/{rat_id}.junc.gz"
    params:
        junc_dir = "{tissue}/splice/junc",
        intermediate = "{tissue}/splice/junc/{rat_id}.junc"
    shell:
        """
        mkdir -p {params.junc_dir}
        regtools junctions extract \
            -a 8 \
            -m 50 \
            -M 500000 \
            -s 0 \
            -o {params.intermediate} \
            {input.bam}
        gzip {params.intermediate}
        """


rule exon_table:
    input:
        "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        "ref/Rattus_norvegicus.Rnor_6.0.99.genes.exons.tsv"
    shell:
        "Rscript src/splice/exon_table.R {input} {output}"


rule splice_bed:
    """Note: fails if intermediate files are already present."""
    input:
        junc = lambda w: expand("{{tissue}}/splice/junc/{rat_id}.junc.gz", rat_id=ids(w.tissue)),
        exons = "ref/Rattus_norvegicus.Rnor_6.0.99.genes.exons.tsv",
        gtf = "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    output:
        bed = "{tissue}/{splice}/{tissue}.leafcutter.bed.gz",
        bedi = "{tissue}/{splice}/{tissue}.leafcutter.bed.gz.tbi",
        clust = "{tissue}/{splice}/{tissue}.leafcutter.clusters_to_genes.txt",
        pcs = "{tissue}/{splice}/{tissue}.leafcutter.PCs.txt",
        groups = "{tissue}/{splice}/{tissue}.leafcutter.phenotype_groups.txt",
    params:
        prefix = "{tissue}",
        script_dir = "src/splice",
        tmpdir = "{tissue}/splice/clust",
    shell:
        """
        mkdir -p {params.tmpdir}
        echo {input.junc} | tr ' ' '\n' > {params.tmpdir}/juncfiles.txt
        python3 src/splice/cluster_prepare_fastqtl.py \
            {params.tmpdir}/juncfiles.txt \
            {input.exons} \
            {input.gtf} \
            {params.prefix} \
            --leafcutter_dir {params.script_dir} \
            --output_dir {params.tmpdir}
        mv -i {params.tmpdir}/{wildcards.tissue}.leafcutter* {params.tmpdir}/..
        rm -r {params.tmpdir}
        """


rule splice_covariates:
    """Compute genotype and splicing PCs and combine."""
    input:
        vcf = "{tissue}/covar/geno.vcf.gz",
        bed = "{tissue}/splice/{tissue}.leafcutter.bed.gz",
    output:
        "{tissue}/splice/{tissue}.covar_splice.txt"
    params:
        n_geno_pcs = 5,
        n_expr_pcs = 10
    shell:
        "Rscript src/covariates.R {input.vcf} {input.bed} {params.n_geno_pcs} {params.n_expr_pcs} {output}"


rule tensorqtl_perm_splice:
    input:
        geno = lambda w: multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{tissue}/splice/{tissue}.covar_splice.txt",
        groups = "{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
    output:
        "{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz"
    params:
        geno_prefix = "{tissue}/geno",
        outdir = "{tissue}/splice"
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
            --groups {input.groups} \
            --mode cis
        """


rule tensorqtl_independent_splice:
    input:
        geno = lambda w: multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{tissue}/splice/{tissue}.covar_splice.txt",
        groups = "{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
        cis = "{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz",
    output:
        "{tissue}/splice/{tissue}_splice.cis_independent_qtl.txt.gz"
    params:
        geno_prefix = "{tissue}/geno",
        outdir = "{tissue}/splice"
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
            --groups {input.groups} \
            --mode cis_independent
        """
