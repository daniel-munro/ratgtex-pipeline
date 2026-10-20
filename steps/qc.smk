localrules:
    qc_sex_concordance,
    qc_star_stats,
    qc_sample_quality_report,
    qc_identity_report,


rule qc_sample_quality_report:
    input:
        stats=expand('{v}/{dataset}/qc/{dataset}.star_stats.tsv', v=VERSION, dataset=DATASETS),
        samples=expand('{v}/{dataset}/rat_ids.txt', v=VERSION, dataset=DATASETS),
    output:
        metrics='reports/sample_quality/sample_metrics.tsv',
        cutoffs='reports/sample_quality/candidate_cutoff_counts.tsv',
        candidates='reports/sample_quality/samples_below_2million.tsv',
        depth_png='reports/sample_quality/uniquely_mapped_depth.png',
        depth_pdf='reports/sample_quality/uniquely_mapped_depth.pdf',
        percent_png='reports/sample_quality/unique_mapping_percent.png',
        percent_pdf='reports/sample_quality/unique_mapping_percent.pdf',
    script:
        '../reports/sample_quality/review_sample_quality.py'


rule qc_identity_report:
    input:
        expand('{v}/{dataset}/qc/rna_to_geno_summary.tsv', v=VERSION, dataset=DATASETS),
    output:
        'reports/identity_qc/issues.tsv',
    script:
        '../reports/identity_qc/summarize.py'


rule qc_mixups_test_snps_vcf:
    """Get subset genotypes for mixup testing.
    Subset to only exon SNPs, SNPs with MAF >= 0.2, and SNPs with <10% missing values.
    """
    input:
        vcf = lambda w: f"geno/{GENOTYPE_PANELS[w.tissue]}.vcf.gz",
        vcfi = lambda w: f"geno/{GENOTYPE_PANELS[w.tissue]}.vcf.gz.tbi",
        regions = "ref/exon_regions.tsv.gz",
    output:
        vcf = "{version}/{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{version}/{tissue}/qc/test_snps.vcf.gz.tbi",
    params:
        out_dir = "{version}/{tissue}/qc"
    shell:
        """
        mkdir -p {params.out_dir}
        bcftools view {input.vcf} \
            --regions-file {input.regions} \
            -Ou | bcftools view \
            --min-af 0.2:minor \
            -i 'F_MISSING<0.1' \
            -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule qc_mixups_ASEReadCounter:
    """Count reads with REF vs. ALT allele in a sample for each SNP."""
    input:
        ref = f"{GENOME_PREFIX}.fa",
        refi = f"{GENOME_PREFIX}.fa.fai",
        refd = f"{GENOME_PREFIX}.dict",
        vcf = "{version}/{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{version}/{tissue}/qc/test_snps.vcf.gz.tbi",
        bam = "{version}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
        bai = "{version}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam.bai",
    output:
        "{version}/{tissue}/qc/test_snps/{rat_id}.readcounts.txt"
    params:
        out_dir = "{version}/{tissue}/qc/test_snps"
    resources:
        mem_mb = 16000,
        runtime = '16h'
    shell:
        """
        mkdir -p {params.out_dir}
        gatk ASEReadCounter \
            -R {input.ref} \
            -I {input.bam} \
            -V {input.vcf} \
            -L {input.vcf} \
            -O {output} \
            --min-depth-of-non-filtered-base 10 \
            --min-mapping-quality 250 \
            --min-base-quality 15
        """


rule qc_mixups_compare_rna_to_vcf:
    """Compare similarity of genotypes and allele counts from RNA-Seq to identify mixups.
    Produces a matrix of RNA-Seq samples x VCF individuals giving similarity
    across the test SNPs. Similarity is mean(1 - abs(RNA_frac - geno_frac)), where RNA_frac
    is fraction of reads with ALT allele for each SNP, and geno_frac is fraction of DNA strands
    with ALT allele (i.e. 0, 0.5, or 1) for each SNP. Also produces a summary table with top
    similarity per RNA-Seq sample to quickly check for mismatches.
    """
    input:
        vcf = "{version}/{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{version}/{tissue}/qc/test_snps.vcf.gz.tbi",
        counts = lambda w: expand("{{version}}/{{tissue}}/qc/test_snps/{rat_id}.readcounts.txt", rat_id=ids(w.tissue)),
        samples = "{version}/{tissue}/rat_ids.txt",
    output:
        matrix = "{version}/{tissue}/qc/rna_to_geno_similarity.tsv",
        summary = "{version}/{tissue}/qc/rna_to_geno_summary.tsv",
    params:
        count_dir = "{version}/{tissue}/qc/test_snps"
    resources:
        runtime = '4h',
        mem_mb = 32000
    shell:
        """
        python3 scripts/qc/rna_to_geno_similarity.py \
            --vcf {input.vcf} \
            --count-dir {params.count_dir} \
            --samples {input.samples} \
            --out-matrix {output.matrix} \
            --out-summary {output.summary}
        """


rule qc_sex_concordance:
    """Compare sex predicted from chrY gene expression to metadata labels."""
    input:
        expr = "{version}/{tissue}/{tissue}.expr.tpm.bed.gz",
        meta = "geno/genotyping_log.csv",
    output:
        "{version}/{tissue}/qc/{tissue}.sex_concordance.txt",
    shell:
        "python3 scripts/qc/sex_concordance.py {input.expr} {input.meta} {output}"


rule qc_star_stats:
    """Parse STAR log files to get stats
    
    The .bai files are not used, but since this rule is used to conveneiently
    generate BAM files, this will produce their .bai files too.
    """
    input:
        bais = lambda w: expand("{{version}}/{{tissue}}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam.bai", rat_id=ids(w.tissue)),
        logs = lambda w: expand("{{version}}/{{tissue}}/star_out/{rat_id}.Log.final.out", rat_id=ids(w.tissue)),
        samples = "{version}/{tissue}/rat_ids.txt",
    output:
        tsv = "{version}/{tissue}/qc/{tissue}.star_stats.tsv",
    params:
        star_dir = "{version}/{tissue}/star_out",
    shell:
        "python3 scripts/qc/stats_from_star_logs.py {params.star_dir} {input.samples} {output}"
