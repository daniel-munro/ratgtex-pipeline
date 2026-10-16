"""Shared round 11.3 panel, tissue hardcalls, and VCFs for alignment/QC/aFC."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(workflow.basedir) / 'scripts'))
from workflow_config import datasets

G = config['genotypes']
GP = f"geno/{G['panel']}"
GI = GP + '.work'
GC = [str(c) for c in G['chromosomes']]
if not GC or [int(c) for c in GC] != sorted(set(int(c) for c in GC)) or any(int(c) not in range(1, 21) for c in GC):
    raise ValueError('Genotype chromosomes must be unique autosomes 1-20 in numeric order')
GSOURCE = G['source_prefix']
GEXONS = GP + '.all_rats_exons'

rule genotypes:
    input:
        expand(GP + '.{ext}', ext=['pgen', 'pvar', 'psam']),
        GP + '.vcf.gz', GP + '.vcf.gz.tbi', GP + '.alleles.tsv.gz',
        expand(config['version'] + '/{tissue}/geno.{ext}', tissue=list(config['tissues']), ext=['bed', 'bim', 'fam'])

rule genotype_sample_map:
    input:
        samples=expand(config['version'] + '/{dataset}/rat_ids.txt', dataset=datasets(config)),
        psams=[GSOURCE.format(chrom=c) + '.psam' for c in GC]
    output:
        keep=GI + '/samples.keep.tsv',
        rename=GI + '/samples.rename.tsv'
    log:
        coverage=GI + '/sample_coverage.tsv',
        counts=GI + '/sample_counts.tsv'
    params:
        datasets=datasets(config)
    script:
        '../scripts/setup/genotype_samples.py'

rule genotype_variant_plan:
    input:
        pvar=lambda w: GSOURCE.format(chrom=w.chrom) + '.pvar',
        fasta=config['ref_genome'],
        fai=config['ref_genome'] + '.fai'
    output:
        extract=GI + '/chr{chrom}/variants.extract.txt',
        reference=GI + '/chr{chrom}/variants.reference.tsv',
        audit=GI + '/chr{chrom}/variants.audit.tsv.gz',
        summary=GI + '/chr{chrom}/variants.summary.json'
    wildcard_constraints:
        chrom='|'.join(GC)
    resources:
        mem_mb=8000
    script:
        '../scripts/setup/genotype_variants.py'

rule genotype_subset_chromosome:
    input:
        source=lambda w: [GSOURCE.format(chrom=w.chrom) + '.' + ext for ext in ['pgen', 'pvar', 'psam']],
        keep=rules.genotype_sample_map.output.keep,
        rename=rules.genotype_sample_map.output.rename,
        extract=rules.genotype_variant_plan.output.extract,
        reference=rules.genotype_variant_plan.output.reference
    output:
        multiext(GI + '/chr{chrom}/ratgtex', '.pgen', '.pvar', '.psam')
    params:
        source=lambda w: GSOURCE.format(chrom=w.chrom),
        prefix=GI + '/chr{chrom}/ratgtex'
    threads: 4
    resources:
        mem_mb=16000,
        tmpdir=GI + '/tmp',
        runtime='4h'
    shell:
        """
        plink2 --pfile {params.source:q} --keep {input.keep:q} \
            --extract {input.extract:q} --ref-allele force {input.reference:q} 2 1 \
            --make-pgen pvar-cols=xheader,maybequal,maybefilter --threads {threads} \
            --memory {resources.mem_mb} --out {params.prefix:q}.selected
        plink2 --pfile {params.prefix:q}.selected --update-ids {input.rename:q} --set-all-var-ids '@:#' \
            --output-chr chrM --sort-vars --make-pgen --threads {threads} \
            --memory {resources.mem_mb} --out {params.prefix:q}
        rm {params.prefix:q}.selected.pgen {params.prefix:q}.selected.pvar {params.prefix:q}.selected.psam
        """

rule genotype_merge_list:
    input:
        expand(GI + '/chr{chrom}/ratgtex.{ext}', chrom=GC, ext=['pgen', 'pvar', 'psam'])
    output:
        GI + '/merge_list.txt'
    run:
        Path(output[0]).write_text(''.join(GI + f'/chr{c}/ratgtex\n' for c in GC))

rule genotype_merge:
    input:
        files=rules.genotype_merge_list.input,
        merge_list=rules.genotype_merge_list.output
    output:
        multiext(GP, '.pgen', '.pvar', '.psam')
    threads: 4
    resources:
        mem_mb=24000,
        tmpdir=GI + '/tmp',
        runtime='8h'
    shell:
        """
        plink2 --pmerge-list {input.merge_list:q} pfile --make-pgen --output-chr chrM \
            --threads {threads} --memory {resources.mem_mb} --out {GP:q}
        """

rule genotype_panel_vcf:
    input:
        multiext(GP, '.pgen', '.pvar', '.psam')
    output:
        vcf=GP + '.vcf.gz',
        index=GP + '.vcf.gz.tbi'
    params:
        prefix=GI + '/panel_export'
    threads: 4
    resources:
        mem_mb=16000,
        tmpdir=GI + '/tmp',
        runtime='8h'
    shell:
        """
        plink2 --pfile {GP:q} --export vcf bgz id-paste=iid --output-chr chrM \
            --threads {threads} --memory {resources.mem_mb} --out {params.prefix:q}
        mv {params.prefix:q}.vcf.gz {output.vcf:q}
        tabix -p vcf {output.vcf:q}
        """

rule genotype_alleles:
    input:
        GP + '.pvar'
    output:
        GP + '.alleles.tsv.gz'
    run:
        import gzip
        sys.path.insert(0, str(Path(workflow.basedir) / 'scripts/setup'))
        from genotype_variants import pvar_rows
        with gzip.open(output[0], 'wt') as handle:
            for chrom, pos, variant, ref, alt in pvar_rows(input[0]):
                handle.write(f'{variant}\t{ref}\t{alt}\n')

rule pgen_to_plink:
    """Apply frequency/count filters within each tissue, then write hardcalls."""
    input:
        panel=multiext(GP, '.pgen', '.pvar', '.psam'),
        samples='{version}/{tissue}/rat_ids.txt'
    output:
        multiext('{version}/{tissue}/geno', '.bed', '.bim', '.fam')
    params:
        prefix='{version}/{tissue}/geno'
    threads: 4
    resources:
        mem_mb=8000,
        tmpdir=GI + '/tmp'
    shell:
        """
        plink2 --pfile {GP:q} --keep {input.samples:q} --make-bed --output-chr chrM \
            --threads {threads} --memory {resources.mem_mb} --out {params.prefix:q}.hardcalls
        plink2 --bfile {params.prefix:q}.hardcalls --maf 0.01 --mac 2 --max-alleles 2 \
            --make-bed --output-chr chrM --threads {threads} --memory {resources.mem_mb} \
            --out {params.prefix:q}
        rm {params.prefix:q}.hardcalls.bed {params.prefix:q}.hardcalls.bim {params.prefix:q}.hardcalls.fam
        """

rule genotype_all_rats_exon_chromosome:
    """Retain every original genotype sample for unexpected-identity searches."""
    input:
        source=lambda w: [GSOURCE.format(chrom=w.chrom) + '.' + ext for ext in ['pgen', 'pvar', 'psam']],
        extract=rules.genotype_variant_plan.output.extract,
        reference=rules.genotype_variant_plan.output.reference,
        regions='ref/exon_regions.tsv.gz'
    output:
        vcf=GI + '/chr{chrom}/all_rats_exons.vcf.gz',
        index=GI + '/chr{chrom}/all_rats_exons.vcf.gz.tbi'
    params:
        source=lambda w: GSOURCE.format(chrom=w.chrom),
        prefix=GI + '/chr{chrom}/all_rats_exons'
    threads: 4
    resources:
        mem_mb=16000,
        tmpdir=GI + '/tmp',
        runtime='4h'
    shell:
        """
        plink2 --pfile {params.source:q} --extract {input.extract:q} \
            --extract-intersect bed1 {input.regions:q} --force-intersect \
            --ref-allele force {input.reference:q} 2 1 --make-pgen pvar-cols=xheader,maybequal,maybefilter \
            --threads {threads} --memory {resources.mem_mb} --out {params.prefix:q}.selected
        plink2 --pfile {params.prefix:q}.selected --set-all-var-ids '@:#' \
            --export vcf bgz id-paste=iid --output-chr chrM \
            --threads {threads} --memory {resources.mem_mb} --out {params.prefix:q}
        rm {params.prefix:q}.selected.pgen {params.prefix:q}.selected.pvar {params.prefix:q}.selected.psam
        tabix -p vcf {output.vcf:q}
        """

rule genotype_all_rats_exons:
    input:
        vcfs=expand(GI + '/chr{chrom}/all_rats_exons.vcf.gz', chrom=GC),
        indexes=expand(GI + '/chr{chrom}/all_rats_exons.vcf.gz.tbi', chrom=GC)
    output:
        vcf=GEXONS + '.vcf.gz',
        index=GEXONS + '.vcf.gz.tbi',
        snps=GEXONS + '.snps.txt'
    threads: 4
    resources:
        mem_mb=8000,
        tmpdir=GI + '/tmp',
        runtime='4h'
    shell:
        """
        bcftools concat {input.vcfs:q} --threads {threads} -Oz -o {output.vcf:q}
        tabix -p vcf {output.vcf:q}
        bcftools query -f '%ID\n' {output.vcf:q} > {output.snps:q}
        """

localrules: collapse_annotation, qc_mixups_exon_regions

rule collapse_annotation:
    """Combine all isoforms of a gene into a single transcript.
    See https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
    """
    input:
        config["ref_anno"]
    output:
        config["ref_anno"].removesuffix(".gtf") + ".genes.gtf"
    shell:
        "python3 scripts/setup/collapse_annotation.py {input} {output}"


rule qc_mixups_exon_regions:
    """Get all exon regions from gene annotations"""
    input:
        config["ref_anno"].removesuffix(".gtf") + ".genes.gtf"
    output:
        "ref/exon_regions.tsv.gz"
    shell:
        """grep -v '^#' {input} | awk '$3=="exon"' | cut -f1,4,5 | gzip -c > {output}"""


