# # Samples for alternative method (test_snps/)
# all_sample_ids = [id for ids in sample_ids.values() for id in ids]
# # all_sample_ids = all_sample_ids[:10]

localrules:
    qc_mixups_exon_regions,
    qc_mixups_test_snps_vcf,
    qc_mixups_compare_rna_to_vcf,

# rule all:
#     input:
#         # expand("readcounts/{sample_id}.readcounts.txt", sample_id=sample_ids["IL"]),
#         # expand("readcounts_mixed/{sample_id}.mixed.readcounts.txt", sample_id=rand_samples),
#         # expand("{group}.counts.txt.gz", group=["Acbc", "IL", "LHB", "PL", "VoLo", "mixed", "trymatch"])
#         "test_snps.counts.txt.gz"


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
    Subset to only exon SNPs, only individuals for the given tissue, SNPs with
    MAF >= 0.2, and SNPs with <10% missing values.
    """
    input:
        vcf = "geno/ratgtex.vcf.gz",
        vcfi = "geno/ratgtex.vcf.gz.tbi",
        samples = "{tissue}/rat_ids.txt",
        regions = "ref/exon_regions.tsv.gz",
    output:
        "{tissue}/qc/test_snps.vcf.gz"
    shell:
        """
        mkdir -p {wildcards.tissue}/qc
        bcftools view {input.vcf} \
            --regions-file {input.regions} \
            --samples-file {input.samples} \
            -Ou | bcftools view \
            --min-af 0.2:minor \
            -i 'F_MISSING<0.1' \
            -Oz -o {output}
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
    output:
        matrix = "{tissue}/qc/rna_to_geno_similarity.tsv",
        summary = "{tissue}/qc/rna_to_geno_summary.tsv",
    params:
        count_dir = "{tissue}/qc/test_snps"
    shell:
        "python3 src/qc_rna_to_geno_similarity.py {input.vcf} {params.count_dir} {output.matrix} {output.summary}"

