rule rsem_index:
    """Generate the index for RSEM."""
    input:
        fasta = f"{GENOME_PREFIX}.fa",
        gtf = f"{ANNO_PREFIX}.gtf"
    output:
        "ref_{rn}/rsem_reference/rsem_reference.transcripts.fa"
    params:
        ref_prefix = "ref_{rn}/rsem_reference/rsem_reference",
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
        ref = "ref_{rn}/rsem_reference/rsem_reference.transcripts.fa",
        bam = "{rn}/{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam"
    output:
        "{rn}/{tissue}/rsem_out/{rat_id}.genes.results.gz"
    params:
        ref_prefix = "ref_{rn}/rsem_reference/rsem_reference",
        out_prefix = "{rn}/{tissue}/rsem_out/{rat_id}",
        paired_end_flag = lambda w: "--paired-end" if config[w.tissue]["paired_end"] else "",
    threads: 16
    resources:
        walltime = 8
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
        rm {params.out_prefix}.isoforms.results
        rm -r {params.out_prefix}.stat
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
        "python3 scripts/collapse_annotation.py {input} {output}"


rule assemble_expression:
    """Combine RSEM output from all samples into log-count and TPM expression tables.
    Also computes inverse-quantile normalized values to use for eQTL mapping. The filtered
    version is to avoid a tensorQTL error on phenotypes with 1 nonzero value.
    """
    input:
        rsem = lambda w: expand("{{rn}}/{{tissue}}/rsem_out/{rat_id}.genes.results.gz", rat_id=ids(w.tissue)),
        anno = f"{ANNO_PREFIX}.genes.gtf"
    output:
        multiext("{rn}/{tissue}/{tissue}.expr.log2.bed", ".gz", ".gz.tbi"),
        multiext("{rn}/{tissue}/{tissue}.expr.tpm.bed", ".gz", ".gz.tbi"),
        multiext("{rn}/{tissue}/{tissue}.expr.iqn.bed", ".gz", ".gz.tbi"),
        multiext("{rn}/{tissue}/{tissue}.expr.iqn.filtered.bed", ".gz", ".gz.tbi")
    params:
        rsem_dir = "{rn}/{tissue}/rsem_out",
        prefix = "{rn}/{tissue}/{tissue}.expr"
    shell:
        """
        python3 scripts/assemble_expression.py {params.rsem_dir} {input.anno} {params.prefix}
        bgzip {params.prefix}.log2.bed
        bgzip {params.prefix}.tpm.bed
        bgzip {params.prefix}.iqn.bed
        bgzip {params.prefix}.iqn.filtered.bed
        tabix {params.prefix}.log2.bed.gz
        tabix {params.prefix}.tpm.bed.gz
        tabix {params.prefix}.iqn.bed.gz
        tabix {params.prefix}.iqn.filtered.bed.gz
        """
