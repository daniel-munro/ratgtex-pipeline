localrules:
    exon_table,
    combine_splice_covariates,


# rat_to_sample = {}
# for tissue in tissues:
#     rat_to_sample[tissue] = {}
#     for sample in sample_ids[tissue]:
#         rat = sample.split('_')[0]
#         rat_to_sample[tissue][rat] = sample
# rat_ids = {}
# for tissue in tissues:
#     rat_ids[tissue] = list(rat_to_sample[tissue].keys())


# def geno_prefix_splice(wildcards):
#     """Genotypes are expression dataset-specific to remove monomorphic SNPs."""
#     tissue = dict(IL="IL", LHb="LHB", NAcc="Acbc", OFC="VoLo", PL="PL")[wildcards.tissue]
#     return f"data/tensorqtl/genotypes/ensembl-gene_inv-quant_{tissue}"


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
        regtools junctions extract -a 8 -m 50 -M 500000 -s 0 {input.bam} -o {params.intermediate}
        gzip {params.intermediate}
        """


rule exon_table:
    input:
        "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        "ref/Rattus_norvegicus.Rnor_6.0.99.genes.exons.tsv"
    shell:
        "Rscript src/splice/exon_table.R {input} {output}"


rule splicing_bed:
    """Note: fails if intermediate files are already present."""
    input:
        junc = lambda w: expand("{{tissue}}/splice/junc/{rat_id}.junc.gz", rat_id=ids[w.tissue]),
        exons = "ref/Rattus_norvegicus.Rnor_6.0.99.genes.exons.tsv",
        gtf = "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf",
    output:
        bed = "{tissue}/{splice}/{tissue}.leafcutter.bed.gz",
        bedi = "{tissue}/{splice}/{tissue}.leafcutter.bed.gz.tbi",
        clust = "{tissue}/{splice}/{tissue}.leafcutter.clusters_to_genes.txt",
        pcs = "{tissue}/{splice}/{tissue}.leafcutter.PCs.txt",
        groups = "{tissue}/{splice}/{tissue}.leafcutter.phenotype_groups.txt",
    params:
        # file_list = "data/splice/{tissue}/juncfiles.txt",
        junc_list = lambda w: expand("../junc/{rat_id}.junc.gz", rat_id=rat_ids[w.tissue]),
        outdir = "{tissue}/splice/clust",
        path_from_outdir = "../../..",
    conda:
        "../envs/splice.yaml"
    shell:
        # paths must be relative to output_dir
        # echo {input.junc} | tr ' ' '\n' > {params.file_list}
        # {params.path_from_outdir}/{params.file_list} \
        # --output_dir data/splice/{wildcards.tissue}
        """
        mkdir -p {params.outdir}
        echo {params.junc_list} | tr ' ' '\n' > {params.outdir}/juncfiles.txt
        cd {params.outdir}
        python2 {params.path_from_outdir}/src/splice/cluster_prepare_fastqtl.py \
            juncfiles.txt \
            {params.path_from_outdir}/{input.exons} \
            {params.path_from_outdir}/{input.gtf} \
            {wildcards.tissue} \
            --leafcutter_dir {params.path_from_outdir}/tools/leafcutter
        mv -i {wildcards.tissue}.leafcutter* ..
        """


rule splicing_covariates:
    """Compute genotype and splicing PCs and combine."""
    input:
        vcf = "{tissue}/covar/geno.vcf.gz",
        bed = "{tissue}/{splice}/{tissue}.leafcutter.bed.gz",
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
