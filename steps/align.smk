localrules:
    # individual_vcf,

rule star_index:
    """Generate the index for STAR.
    A different index is generated for each read length.
    """
    input:
        fasta = f"{GENOME_PREFIX}.fa",
        gtf = f"{ANNO_PREFIX}.gtf",
    output:
        # Among others:
        "ref/star_index_{read_length}/SAindex"
    params:
        outdir = "ref/star_index_{read_length}",
        overhang = lambda w: int(w.read_length) - 1
    threads: 8
    resources:
        mem_mb = 60000,
        runtime = '4h'
    shell:
        """
        mkdir -p {params.outdir}
        STAR --runMode genomeGenerate \
            --outTmpDir {params.outdir}/tmp \
            --outFileNamePrefix {params.outdir}/STAR_ \
            --genomeDir {params.outdir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --runThreadN {threads}
        """


rule individual_vcf:
    """Get an individual-specific VCF file.
    This is used by STAR to consider the individual's variants for better alignment.
    """
    input:
        "geno/{geno_dataset}.vcf.gz"
    output:
        "geno/individual/{geno_dataset}/{rat_id}.vcf.gz"
    params:
        outdir = "geno/individual/{geno_dataset}"
    shell:
        """
        mkdir -p {params.outdir}
        bcftools view -s {wildcards.rat_id} --min-ac=1 -O z -o {output} {input}
        """


def fastq_input(wildcards) -> list:
    """Get the list of FASTQ file paths for a sample.
    
    Returns a list of paths. For paired-end samples, file pairs are adjacent
    in the list, e.g. [runA_1.fq.gz runA_2.fq.gz runB_1.fq.gz runB_2.fq.gz].
    For single-end samples, returns just the list of files.
    """
    paths, is_paired = FASTQ_MAPS[wildcards.tissue][wildcards.rat_id]
    if is_paired:
        result = []
        for i in range(len(paths[0])):
            result.append(paths[0][i])
            result.append(paths[1][i])
        return result
    else:
        return paths


def fastq_star_param(wildcards):
    """Get a string listing the fastq input files, with two lists if paired_end.
    
    This is the string supplied directly to the STAR command.
    """
    paths, is_paired = FASTQ_MAPS[wildcards.tissue][wildcards.rat_id]
    if is_paired:
        return ' '.join([','.join(paths[0]), ','.join(paths[1])])
    else:
        return ','.join(paths)


def read_groups(wildcards, input):
    """Include read group in BAM for MarkDuplicates (though we currently aren't using that)."""
    _, is_paired = FASTQ_MAPS[wildcards.tissue][wildcards.rat_id]
    n_fastqs = len(input.fastq) // 2 if is_paired else len(input.fastq)
    rgs = expand("ID:{fq} SM:{sam}", fq=range(n_fastqs), sam=wildcards.rat_id)
    return " , ".join(rgs)


rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        fastq = fastq_input,
        vcf = lambda w: f"geno/individual/{config['tissues'][w.tissue]['geno_dataset']}/{w.rat_id}.vcf.gz",
        index = lambda w: f"ref/star_index_{config['tissues'][w.tissue]['read_length']}/SAindex"
    output:
        # RSEM requires transcriptome-sorted BAM.
        coord = "{version}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bam = "{version}/{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam",
        log = "{version}/{tissue}/star_out/{rat_id}.Log.final.out",
    params:
        fastq_list = fastq_star_param,
        index_dir = lambda w: f"ref/star_index_{config['tissues'][w.tissue]['read_length']}",
        prefix = "{version}/{tissue}/star_out/{rat_id}.",
        read_groups = read_groups,
        out_dir = "{version}/{tissue}/star_out"
    threads: 16
    resources:
        mem_mb = 60000,
        runtime = '8h'
    shell:
        """
        mkdir -p {params.out_dir}
        STAR --runMode alignReads \
            --runThreadN {threads} \
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


rule index_bam:
    """Index a BAM file."""
    input:
        "{version}/{tissue}/star_out/{basename}.bam",
    output:
        "{version}/{tissue}/star_out/{basename}.bam.bai",
    params:
        add_threads = 8 - 1,
    threads: 8
    shell:
        # It expects the number of *additional* threads to use beyond the first.
        """
        samtools index -@ {params.add_threads} {input}
        """


