localrules:
    # individual_vcf,

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
        overhang = read_length - 1
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
        --runThreadN {resources.cpus}
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


def fastqs(tissue: str, rat_id: str, paired: bool) -> list:
    """Get the list of FASTQ file paths for a sample.
    If paired_end is True, return a list of two lists of corresponding paths.
    """
    paths = [[], []] if paired_end else []
    with open(f"{tissue}/fastq_map.txt", "r") as f:
        for line in f.read().splitlines():
            if paired_end:
                fastq1, fastq2, r_id = line.split("\t")
                if r_id == rat_id:
                    paths[0].append(str(fastq_path / fastq1))
                    paths[1].append(str(fastq_path / fastq2))
            else:
                fastq, r_id = line.split("\t")
                if r_id == rat_id:
                    paths.append(str(fastq_path / fastq))
    return paths


def fastq_input(wildcards):
    """Get the FASTQ input file/s for a sample.
    If paired_end is True, concatenate into one list of files.
    """
    files = fastqs(wildcards.tissue, wildcards.rat_id, paired_end)
    return files[0] + files[1] if paired_end else files


def fastq_param(wildcards, input):
    """Get a string listing the fastq input files, with two lists if paired_end."""
    if paired_end:
        files = fastqs(wildcards.tissue, wildcards.rat_id, True)
        return " ".join([",".join(files[0]), ",".join(files[1])])
    else:
        return ",".join(input.fastq)


def read_groups(wildcards, input):
    """Include read group in BAM for MarkDuplicates (though we currently aren't using that)."""
    n_fastqs = len(input.fastq) // 2 if paired_end else len(input.fastq)
    rgs = expand("ID:{fq} SM:{sam}", fq=range(n_fastqs), sam=wildcards.rat_id)
    return " , ".join(rgs)


rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        # fastq = lambda w: list(fastqs.loc[fastqs["rat_id"] == w.rat_id, "path"]),
        fastq = fastq_input,
        vcf = "geno/individual/{rat_id}.vcf.gz",
        index = "{tissue}/star_index/SAindex"
    output:
        # RSEM requires transcriptome-sorted BAM.
        coord = "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bam = "{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam"
    params:
        fastq_list = fastq_param,
        index_dir = "{tissue}/star_index",
        prefix = "{tissue}/star_out/{rat_id}.",
        read_groups = read_groups,
    resources:
        mem_mb = 60000,
        cpus = 16,
        walltime = 8
    shell:
        """
        mkdir -p {wildcards.tissue}/star_out
        STAR --runMode alignReads \
            --runThreadN {resources.cpus} \
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
