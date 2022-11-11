localrules:
    qc_mixups_exon_regions,
    qc_mixups_test_snps_vcf,
    # qc_mixups_compare_rna_to_vcf,
    qc_sex_concordance,


rule qc_mixups_exon_regions:
    """Get all exon regions from gene annotations"""
    input:
        f"{ANNO_PREFIX}.genes.gtf"
    output:
        "ref_{rn}/exon_regions.tsv.gz"
    shell:
        """grep -v '^#' {input} | awk '$3=="exon"' | cut -f1,4,5 | gzip -c > {output}"""


rule qc_mixups_test_snps_vcf:
    """Get subset genotypes for mixup testing.
    Subset to only exon SNPs, SNPs with MAF >= 0.2, and SNPs with <10% missing values.
    """
    input:
        vcf = lambda w: f"geno_{RN}/{config[w.tissue]['geno_dataset']}.vcf.gz",
        vcfi = lambda w: f"geno_{RN}/{config[w.tissue]['geno_dataset']}.vcf.gz.tbi",
        # samples = "{rn}/{tissue}/rat_ids.txt",
        regions = "ref_{rn}/exon_regions.tsv.gz",
    output:
        vcf = "{rn}/{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{rn}/{tissue}/qc/test_snps.vcf.gz.tbi",
    shell:
            # --samples-file {input.samples} \
        """
        mkdir -p {wildcards.rn}/{wildcards.tissue}/qc
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
        vcf = "{rn}/{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{rn}/{tissue}/qc/test_snps.vcf.gz.tbi",
        # bam = "{rn}/{tissue}/markdup_out/{rat_id}.bam",
        bam = "{rn}/{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
    output:
        "{rn}/{tissue}/qc/test_snps/{rat_id}.readcounts.txt"
    resources:
        walltime = 16
    shell:
        """
        mkdir -p {wildcards.rn}/{wildcards.tissue}/qc/test_snps
        gatk ASEReadCounter \
            -R {input.ref} \
            -I {input.bam} \
            -V {input.vcf} \
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
        vcf = "{rn}/{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{rn}/{tissue}/qc/test_snps.vcf.gz.tbi",
        counts = lambda w: expand("{{rn}}/{{tissue}}/qc/test_snps/{rat_id}.readcounts.txt", rat_id=ids(w.tissue)),
        samples = "{rn}/{tissue}/rat_ids.txt",
    output:
        matrix = "{rn}/{tissue}/qc/rna_to_geno_similarity.tsv",
        summary = "{rn}/{tissue}/qc/rna_to_geno_summary.tsv",
    params:
        count_dir = "{rn}/{tissue}/qc/test_snps"
    resources:
        walltime = 4,
        mem_mb = 16000
    shell:
        """
        python3 scripts/qc_rna_to_geno_similarity.py \
            {input.vcf} \
            {params.count_dir} \
            {input.samples} \
            {output.matrix} \
            {output.summary}
        """


# rule qc_mixups_test_snps_vcf_all_rats:
#     """Get subset genotypes for all available rats for additional mixup testing.
#     Subset to the same set of test SNPs that were counted in the RNA samples.
#     """
#     input:
#         all_rat_vcf = "geno/all_rats_exons.vcf.gz",
#         all_rat_vcfi = "geno/all_rats_exons.vcf.gz.tbi",
#         test_snp_vcf = "{tissue}/qc/test_snps.vcf.gz",
#         test_snp_vcfi = "{tissue}/qc/test_snps.vcf.gz",
#     output:
#         "{tissue}/qc/test_snps.all_rats.vcf.gz"
#     shell:
#         """
#         bcftools view {input.all_rat_vcf} \
#             --targets-file {input.test_snp_vcf} \
#             -Oz -o {output}
#         """


rule qc_mixups_compare_to_all_rats:
    """For samples without a genotype match, check for matches in all rats.
    After running qc_mixups_compare_rna_to_vcf, some mismatched RNA samples might still
    not have a genotype match. This rule will compare their test SNPs to those of all
    6000+ rats we have genotypes for to see if any match. It's good to also include at
    least one RNA sample ID that did match as a positive control (as long as it's included
    in the all-rat VCF).
    """
    input:
        vcf = "geno_{rn}/all_rats_exons.vcf.gz",
        vcfi = "geno_{rn}/all_rats_exons.vcf.gz.tbi",
        samples = "{rn}/{tissue}/qc/samples_without_matches.txt",
    output:
        "{rn}/{tissue}/qc/all_rats_summary.tsv"
    params:
        count_dir = "{rn}/{tissue}/qc/test_snps"
    resources:
        walltime = 4,
        mem_mb = 16000
    shell:
        """
        python3 scripts/qc_rna_to_geno_all_rats.py \
            {input.vcf} \
            {params.count_dir} \
            {input.samples} \
            {output}
        """

rule qc_sex_concordance:
    """Compare sex predicted from chrY gene expression to metadata labels."""
    input:
        expr = "{rn}/{tissue}/{tissue}.expr.tpm.bed.gz",
        meta = "geno_{rn}/genotyping_log.csv",
    output:
        "{rn}/{tissue}/qc/{tissue}.sex_concordance.txt",
    shell:
        "python3 scripts/qc_sex_concordance.py {input.expr} {input.meta} {output}"
