localrules:
    collapse_annotation,


rule rsem_index:
    """Generate the index for RSEM."""
    input:
        fasta = f"{GENOME_PREFIX}.fa",
        gtf = f"{ANNO_PREFIX}.gtf"
    output:
        "ref/rsem_reference/rsem_reference.transcripts.fa"
    params:
        ref_prefix = "ref/rsem_reference/rsem_reference",
    threads: 16
    shell:
        """
        rsem-prepare-reference \
            {input.fasta} \
            {params.ref_prefix} \
            --gtf {input.gtf} \
            --num-threads {threads}
        """


rule rsem:
    """Quantify expression from a BAM file."""
    input:
        ref = "ref/rsem_reference/rsem_reference.transcripts.fa",
        bam = "{version}/{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam"
    output:
        "{version}/{tissue}/rsem_out/{rat_id}.genes.results.gz"
    params:
        ref_prefix = "ref/rsem_reference/rsem_reference",
        out_prefix = "{version}/{tissue}/rsem_out/{rat_id}",
        paired_end_flag = lambda w: "--paired-end" if FASTQ_MAPS[w.tissue][w.rat_id][1] else "",
    threads: 16
    resources:
        runtime = '8h'
    shell:
        """
        rsem-calculate-expression \
            {params.paired_end_flag} \
            --num-threads {threads} \
            --quiet \
            --estimate-rspd \
            --no-bam-output \
            --alignments {input.bam} \
            {params.ref_prefix} \
            {params.out_prefix}
        gzip {params.out_prefix}.genes.results
        gzip {params.out_prefix}.isoforms.results
        """


rule collapse_annotation:
    """Combine all isoforms of a gene into a single transcript.
    See https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
    """
    input:
        f"{ANNO_PREFIX}.gtf"
    output:
        f"{ANNO_PREFIX}.genes.gtf"
    shell:
        "python3 scripts/setup/collapse_annotation.py {input} {output}"


rule assemble_expression:
    """Combine RSEM output from all samples into log-count and TPM expression tables.
    Also computes inverse-quantile normalized values to use for eQTL mapping. The filtered
    version is to avoid a tensorQTL error on phenotypes with 1 nonzero value.
    """
    input:
        rsem = lambda w: expand("{{version}}/{{tissue}}/rsem_out/{rat_id}.genes.results.gz", rat_id=ids(w.tissue)),
        samples = "{version}/{tissue}/rat_ids.txt",
        anno = f"{ANNO_PREFIX}.genes.gtf"
    output:
        multiext("{version}/{tissue}/{tissue}.expr.log2.bed", ".gz", ".gz.tbi"),
        multiext("{version}/{tissue}/{tissue}.expr.tpm.bed", ".gz", ".gz.tbi"),
        multiext("{version}/{tissue}/{tissue}.expr.iqn.bed", ".gz", ".gz.tbi"),
        multiext("{version}/{tissue}/{tissue}.expr.iqn.filtered.bed", ".gz", ".gz.tbi")
    params:
        rsem_dir = "{version}/{tissue}/rsem_out",
        prefix = "{version}/{tissue}/{tissue}.expr"
    shell:
        """
        python3 scripts/assemble_expression.py {params.rsem_dir} {input.samples} {input.anno} {params.prefix}
        bgzip {params.prefix}.log2.bed
        bgzip {params.prefix}.tpm.bed
        bgzip {params.prefix}.iqn.bed
        bgzip {params.prefix}.iqn.filtered.bed
        tabix {params.prefix}.log2.bed.gz
        tabix {params.prefix}.tpm.bed.gz
        tabix {params.prefix}.iqn.bed.gz
        tabix {params.prefix}.iqn.filtered.bed.gz
        """
