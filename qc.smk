localrules:
    qc_mixups_exon_regions,
    qc_mixups_test_snps_vcf,
    qc_mixups_compare_rna_to_vcf,


rule qc_mixups_exon_regions:
    """Get all exon regions from gene annotations"""
    input:
        "ref/Rattus_norvegicus.Rnor_6.0.99.genes.gtf"
    output:
        "ref/exon_regions.tsv.gz"
    shell:
        """grep -v '^#' {input} | awk '$3=="exon"' | cut -f1,4,5 | gzip -c > {output}"""


rule qc_mixups_test_snps_vcf:
    """Get subset genotypes for mixup testing.
    Subset to only exon SNPs, SNPs with MAF >= 0.2, and SNPs with <10% missing values.
    """
    input:
        vcf = f"geno/{geno_dataset}.vcf.gz",
        vcfi = f"geno/{geno_dataset}.vcf.gz.tbi",
        # samples = "{tissue}/rat_ids.txt",
        regions = "ref/exon_regions.tsv.gz",
    output:
        vcf = "{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{tissue}/qc/test_snps.vcf.gz.tbi",
    shell:
            # --samples-file {input.samples} \
        """
        mkdir -p {wildcards.tissue}/qc
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
        ref = "ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
        refi = "ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.fai",
        refd = "ref/Rattus_norvegicus.Rnor_6.0.dna.toplevel.dict",
        vcf = "{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{tissue}/qc/test_snps.vcf.gz.tbi",
        # bam = "{tissue}/markdup_out/{rat_id}.bam",
        bam = "{tissue}/star_out/{rat_id}.Aligned.sortedByCoord.out.bam",
    output:
        "{tissue}/qc/test_snps/{rat_id}.readcounts.txt"
    resources:
        walltime = 16
    shell:
        """
        mkdir -p {wildcards.tissue}/qc/test_snps
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
        vcf = "{tissue}/qc/test_snps.vcf.gz",
        vcfi = "{tissue}/qc/test_snps.vcf.gz.tbi",
        counts = lambda w: expand("{{tissue}}/qc/test_snps/{rat_id}.readcounts.txt", rat_id=ids(w.tissue)),
        samples = "{tissue}/rat_ids.txt",
    output:
        matrix = "{tissue}/qc/rna_to_geno_similarity.tsv",
        summary = "{tissue}/qc/rna_to_geno_summary.tsv",
    params:
        count_dir = "{tissue}/qc/test_snps"
    shell:
        "python3 src/qc_rna_to_geno_similarity.py {input.vcf} {params.count_dir} {input.samples} {output.matrix} {output.summary}"


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
        vcf = "geno/all_rats_exons.vcf.gz",
        vcfi = "geno/all_rats_exons.vcf.gz.tbi",
        samples = "{tissue}/qc/samples_without_matches.txt",
    output:
        "{tissue}/qc/all_rats_summary.tsv"
    params:
        count_dir = "{tissue}/qc/test_snps"
    resources:
        walltime = 4,
        mem_mb = 16000
    shell:
        "python3 src/qc_rna_to_geno_all_rats.py {input.vcf} {params.count_dir} {input.samples} {output}"
