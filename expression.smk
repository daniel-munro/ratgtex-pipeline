rule rsem_index:
    """Generate the index for RSEM."""
    input:
        fasta = "ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        gtf = "ref/Rattus_norvegicus.Rnor_6.0.99.gtf"
    output:
        "ref/rsem_reference/rsem_reference.transcripts.fa"
    resources:
        cpus = 16
    shell:
        """
        rsem-prepare-reference \
        {input.fasta} \
        ref/rsem_reference/rsem_reference \
        --gtf {input.gtf} \
        --num-threads {resources.cpus}
        """


rule rsem:
    """Quantify expression from a BAM file."""
    input:
        ref = "ref/rsem_reference/rsem_reference.transcripts.fa",
        bam = "{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam"
    output:
        "{tissue}/rsem_out/{rat_id}.genes.results.gz"
    params:
        ref_prefix = "ref/rsem_reference/rsem_reference",
        out_prefix = "{tissue}/rsem_out/{rat_id}",
        paired_end_flag = "--paired-end" if paired_end else "",
    resources:
        cpus = 16,
        walltime = 8
    shell:
        """
        rsem-calculate-expression \
        {params.paired_end_flag} \
        --num-threads {resources.cpus} \
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
        "ref/Rattus_norvegicus.Rnor_6.0.99.gtf"
    output:
        "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    shell:
        "python3 src/collapse_annotation.py {input} {output}"


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
