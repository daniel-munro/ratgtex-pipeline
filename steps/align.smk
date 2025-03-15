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
        "ref_{rn}/star_index_{read_length}/SAindex"
    params:
        outdir = "ref_{rn}/star_index_{read_length}",
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
        "geno_{rn}/{geno_dataset}.vcf.gz"
    output:
        "geno_{rn}/individual/{geno_dataset}/{rat_id}.vcf.gz"
    params:
        outdir = "geno_{rn}/individual/{geno_dataset}"
    shell:
        """
        mkdir -p {params.outdir}
        bcftools view -s {wildcards.rat_id} --min-ac=1 -O z -o {output} {input}
        """


def fastqs(tissue: str, rat_id: str, paired: bool) -> list:
    """Get the list of FASTQ file paths for a sample.
    If paired is True, return a list of two lists of corresponding paths.
    """
    fastq_path = Path(config["tissues"][tissue]["fastq_path"])
    paths = [[], []] if paired else []
    with open(f"{RN}/{tissue}/fastq_map.txt", "r") as f:
        for line in f.read().splitlines():
            if paired:
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
    If paired_end is True for this tissue, concatenate into one list of files.
    """
    paired_end = config["tissues"][wildcards.tissue]["paired_end"]
    files = fastqs(wildcards.tissue, wildcards.rat_id, paired_end)
    return files[0] + files[1] if paired_end else files


def fastq_param(wildcards, input):
    """Get a string listing the fastq input files, with two lists if paired_end."""
    if config["tissues"][wildcards.tissue]["paired_end"]:
        files = fastqs(wildcards.tissue, wildcards.rat_id, True)
        return " ".join([",".join(files[0]), ",".join(files[1])])
    else:
        return ",".join(input.fastq)


def read_groups(wildcards, input):
    """Include read group in BAM for MarkDuplicates (though we currently aren't using that)."""
    paired_end = config["tissues"][wildcards.tissue]["paired_end"]
    n_fastqs = len(input.fastq) // 2 if paired_end else len(input.fastq)
    rgs = expand("ID:{fq} SM:{sam}", fq=range(n_fastqs), sam=wildcards.rat_id)
    return " , ".join(rgs)


rule star_align:
    """Align RNA-Seq reads for a sample using STAR."""
    input:
        fastq = fastq_input,
        vcf = lambda w: f"geno_{RN}/individual/{config['tissues'][w.tissue]['geno_dataset']}/{w.rat_id}.vcf.gz",
        index = lambda w: f"ref_{RN}/star_index_{config['tissues'][w.tissue]['read_length']}/SAindex"
    output:
        # RSEM requires transcriptome-sorted BAM.
        coord = "{rn}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bam = "{rn}/{tissue}/star_out/{rat_id}.Aligned.toTranscriptome.out.bam"
    params:
        fastq_list = fastq_param,
        index_dir = lambda w: f"ref_{RN}/star_index_{config['tissues'][w.tissue]['read_length']}",
        prefix = "{rn}/{tissue}/star_out/{rat_id}.",
        read_groups = read_groups,
    threads: 16
    resources:
        mem_mb = 60000,
        runtime = '20h'
    shell:
        """
        mkdir -p {wildcards.rn}/{wildcards.tissue}/star_out
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
        "{rn}/{tissue}/star_out/{basename}.bam",
    output:
        "{rn}/{tissue}/star_out/{basename}.bam.bai",
    params:
        add_threads = 8 - 1,
    threads: 8
    shell:
        # It expects the number of *additional* threads to use beyond the first.
        """
        samtools index -@ {params.add_threads} {input}
        """


