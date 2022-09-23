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
        "Rscript scripts/splice/exon_table.R {input} {output}"


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
        script_dir = "scripts/splice",
        tmpdir = "{tissue}/splice/clust",
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
        "Rscript scripts/covariates.R {input.vcf} {input.bed} {params.n_geno_pcs} {params.n_expr_pcs} {output}"


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
        python3 scripts/run_tensorqtl.py \
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
        geno = multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
        bed = "{tissue}/splice/{tissue}.leafcutter.bed.gz",
        bedi = "{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
        covar = "{tissue}/splice/{tissue}.covar_splice.txt",
    output:
        "{tissue}/splice/{tissue}_splice.trans_qtl_pairs.txt.gz"
    params:
        geno_prefix = "{tissue}/geno",
        outdir = "{tissue}/splice",
        out_prefix = "{tissue}_splice",
    resources:
        walltime = 12,
        # partition = "--partition=gpu",
    shell:
        # batch_size set due to "RuntimeError: CUDA out of memory"
        # tensorQTL doesn't use group info for trans mapping
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


## This gives an error since the groups files are out of order, so I'll omit for now.
# rule tensorqtl_nominal_splice:
#     """Get summary statistics for all tested cis-window SNPs per gene."""
#     input:
#         geno = multiext("{tissue}/geno", ".bed", ".bim", ".fam"),
#         bed = "{tissue}/splice/{tissue}.leafcutter.bed.gz",
#         bedi = "{tissue}/splice/{tissue}.leafcutter.bed.gz.tbi",
#         covar = "{tissue}/splice/{tissue}.covar_splice.txt",
#         groups = "{tissue}/splice/{tissue}.leafcutter.phenotype_groups.txt",
#     output:
#         expand("{{tissue}}/splice/nominal/{{tissue}}_splice.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21))
#     params:
#         geno_prefix = "{tissue}/geno",
#         outdir = "{tissue}/splice/nominal",
#         out_prefix = "{tissue}_splice",
#     resources:
#         walltime = 12,
#         # partition = "--partition=gpu",
#     shell:
#         # module load cuda
#         """
#         mkdir -p {params.outdir}
#         python3 -m tensorqtl \
#             {params.geno_prefix} \
#             {input.bed} \
#             {params.out_prefix} \
#             --covariates {input.covar} \
#             --phenotype_groups {input.groups} \
#             --output_dir {params.outdir} \
#             --mode cis_nominal
#         """


# rule tensorqtl_all_signif_splice:
#     """Extract all significant cis SNP-gene pairs."""
#     input:
#         perm = "{tissue}/splice/{tissue}_splice.cis_qtl.txt.gz",
#         nom = expand("{{tissue}}/splice/nominal/{{tissue}}_splice.cis_qtl_pairs.{chrn}.parquet", chrn=range(1, 21)),
#     output:
#         "{tissue}/splice/{tissue}_splice.cis_qtl_signif.txt.gz"
#     params:
#         nom_prefix = "{tissue}/splice/nominal/{tissue}_splice",
#     shell:
#         "python3 scripts/tensorqtl_all_signif.py {input.perm} {params.nom_prefix} {output}"

