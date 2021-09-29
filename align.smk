localrules:
    individual_vcf,

rule star_index:
    """Generate the index for STAR.
    The index is tissue-specific since the read length can differ among datasets.
    """
    input:
        fasta = "ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        gtf = "ref/Rattus_norvegicus.Rnor_6.0.99.gtf"
    output:
        # Among others:
        "{tissue}/star_index/SAindex"
    params:
        outdir = "{tissue}/star_index",
        overhang = lambda w: read_length - 1
    # threads: 8
    resources:
        cpus = 8
    shell:
        """
        STAR --runMode genomeGenerate \
        --genomeDir {params.outdir} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang {params.overhang} \
        --runThreadN 8
        """


rule individual_vcf:
    """Get an individual-specific VCF file.
    This is used by STAR to consider the individual's variants for better alignment.
    """
    input:
        "geno/ratgtex.vcf.gz"
    output:
        "geno/individual/{rat_id}.vcf.gz"
    shell:
        "bcftools view -s {wildcards.rat_id} --min-ac=1 -O z -o {output} {input}"


def fastqs(wildcards):
    """Get the list of FASTQ file paths for a sample."""
    paths = []
    with open(f"{wildcards.tissue}/fastq_map.txt", "r") as f:
        for line in f.read().splitlines():
            fastq, rat_id = line.split("\t")
            if rat_id == wildcards.rat_id:
                paths.append(str(fastq_path / fastq))
    return paths


def read_groups(wildcards, input):
    """Include read group in BAM for MarkDuplicates (though we currently aren't using that)."""
    rgs = expand("ID:{fq} SM:{sam}", fq=input.fastq, sam=wildcards.rat_id)
    return " , ".join(rgs)


rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        # fastq = lambda w: list(fastqs.loc[fastqs["rat_id"] == w.rat_id, "path"]),
        fastq = fastqs,
        vcf = "geno/individual/{rat_id}.vcf.gz",
        index = "{tissue}/star_index/SAindex"
    output:
        # RSEM requires transcriptome-sorted BAM.
        coord = "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bam = "{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam"
    params:
        fastq_list = lambda wildcards, input: ",".join(input.fastq),
        index_dir = "{tissue}/star_index",
        prefix = "{tissue}/star_out/{rat_id}.",
        read_groups = read_groups,
    resources:
        mem_mb = 60000,
        cpus = 16,
        walltime = 2
    shell:
        """
        STAR --runMode alignReads \
        --runThreadN 16 \
        --genomeDir {params.index_dir} \
        --readFilesIn {params.fastq_list} \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --varVCFfile <(zcat {input.vcf}) \
        --waspOutputMode SAMtag \
        --quantMode TranscriptomeSAM \
        --outSAMstrandField intronMotif \
        --outSAMattrRGline {params.read_groups} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFileNamePrefix {params.prefix}
        """


# WASP-filtering and MarkDuplicates steps were for allele-specific expression,
# but for now we don't need them for this pipeline.

# rule filter_wasp:
#     input:
#         "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam"
#     output:
#         "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.wasp.bam"
#     resources:
#         cpus = 16,
#         walltime = 6
#     shell:
#         """
#         samtools view -H {input} > {output}.sam
#         samtools view {input} | grep -v "vW:i:[2-7]" >> {output}.sam
#         samtools view -1 {output}.sam --threads 16 > {output}
#         rm {output}.sam
#         """


# rule mark_duplicates:
#     input:
#         "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.wasp.bam"
#     output:
#         bam = "{tissue}/markdup_out/{rat_id}.bam",
#         metrics = "{tissue}/markdup_out/{rat_id}.marked_dup_metrics.txt"
#     shell:
#         # MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=500
#         """
#         picard MarkDuplicates \
#         I={input} \
#         O={output.bam} \
#         M={output.metrics} \
#         ASSUME_SORT_ORDER=coordinate \
#         PROGRAM_RECORD_ID=null \
#         TMP_DIR=$TMPDIR \
#         MAX_RECORDS_IN_RAM=2000000
#         """
