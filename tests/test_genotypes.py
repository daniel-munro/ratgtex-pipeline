"""Genotype identity failures and a real PLINK/bcftools miniature workflow."""
import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pandas as pd
import pysam
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts/setup'))
from genotype_samples import prepare_samples


class GenotypeTests(unittest.TestCase):
    def setUp(self):
        (ROOT / '.snakemake').mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / '.snakemake')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)

    def test_sample_coverage_and_ambiguity(self):
        w = self.work
        sample = w / 'rats.txt'
        sample.write_text('rat_A\nabsent\n')
        psam = w / 'samples.psam'
        psam.write_text('#IID\tSEX\nRattaca1_rat_A_x\t1\n')
        outputs = [w / name for name in ['keep', 'rename', 'coverage', 'counts']]
        with self.assertRaisesRegex(ValueError, 'Missing or ambiguous'):
            prepare_samples([sample], ['dataset'], [psam], *outputs)
        coverage = pd.read_csv(outputs[2], sep='\t')
        self.assertEqual(list(coverage.status), ['present', 'missing'])
        sample.write_text('rat_A\n')
        psam.write_text(psam.read_text() + 'Riptide2_rat_A_x\t1\n')
        with self.assertRaisesRegex(ValueError, 'Missing or ambiguous'):
            prepare_samples([sample], ['dataset'], [psam], *outputs)
        self.assertEqual(pd.read_csv(outputs[2], sep='\t').status[0], 'ambiguous')
        other = w / 'other.psam'
        other.write_text('#IID\tSEX\nRattaca1_rat_A_x\t2\n')
        with self.assertRaisesRegex(ValueError, 'PSAM differs'):
            prepare_samples([sample], ['dataset'], [psam, other], *outputs)

    def test_miniature_workflow(self):
        w = self.work
        for folder in ['input', 'ref', 'vtest/First', 'vtest/Second', 'geno']:
            (w / folder).mkdir(parents=True, exist_ok=True)
        (w / 'vtest/First/rat_ids.txt').write_text('rat_A\nrat_B\n')
        (w / 'vtest/Second/rat_ids.txt').write_text('rat_A\nrat_B\nrat_C\n')
        fasta = w / 'ref/genome.fa'
        fasta.write_text('>chr1\n' + 'A' * 100 + '\n>chr2\n' + 'A' * 100 + '\n')
        pysam.faidx(str(fasta))
        with gzip.open(w / 'ref/exon_regions.tsv.gz', 'wt') as handle:
            handle.write('chr1\t1\t30\nchr2\t1\t30\n')
        ids = ['Rattaca1_rat_A_x', 'ddGBS_HSR_1_rat_B', 'LSWrnaseq_rat_C_LSWNNNNN', 'Riptide1_unexpected_x']
        for chrom in [1, 2]:
            vcf = w / f'input/chr{chrom}.vcf'
            header = ('##fileformat=VCFv4.2\n' + f'##contig=<ID={chrom},length=100>\n'
                      '##INFO=<ID=EAF,Number=1,Type=Float,Description="Original cohort frequency">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '##FORMAT=<ID=DS,Number=1,Type=Float,Description="Dosage">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(ids) + '\n')
            records = []
            for pos, name, ref, alt in [(10, 'keep', 'A', 'G'), (20, 'swap', 'G', 'A'),
                                        (30, 'mismatch', 'C', 'T'), (40, 'dup1', 'A', 'G'),
                                        (40, 'dup2', 'A', 'C'), (50, 'indel', 'A', 'AG'),
                                        (60, 'mono', 'A', 'G')]:
                calls = ['0/0:0', '1/1:2', '0/1:0.8', '0/1:1'] if name != 'mono' else ['0/0:0', '0/0:0', '1/1:2', '0/1:1']
                records.append(f'{chrom}\t{pos}\t{name}{chrom}\t{ref}\t{alt}\t.\tPASS\tEAF=0.7\tGT:DS\t' + '\t'.join(calls) + '\n')
            vcf.write_text(header + ''.join(records))
            self.command(['plink2', '--vcf', str(vcf), 'dosage=DS', '--make-pgen', '--out', str(w / f'input/chr{chrom}')])
        config = dict(version='vtest', genotypes=dict(panel='test', source_prefix='input/chr{chrom}', chromosomes=[1, 2]),
                      ref_genome='ref/genome.fa', ref_anno='ref/annotation.gtf',
                      tissues={'First': ['First'], 'Second': ['Second']},
                      datasets={name: dict(read_length=100, fastq_path='reads') for name in ['First', 'Second']},
                      modalities={'expression': dict(grouped=False, producer='expression')},
                      cross_modality=False, default_stage='align')
        config_path = w / 'config.yaml'
        config_path.write_text(yaml.safe_dump(config))
        targets = ['genotypes', 'geno/test.all_rats_exons.vcf.gz', 'geno/test.all_rats_exons.snps.txt',
                   'vtest/First/geno.bed', 'vtest/Second/geno.bed']
        # Genotypes use the main workflow even before any FASTQ maps exist.
        self.command(['snakemake', '-s', str(ROOT / 'Snakefile'), '--configfile', str(config_path),
                      '--profile', 'slurm', '--jobs', '2', '--dry-run', '--', 'genotypes'])
        self.command(['snakemake', '-s', str(ROOT / 'Snakefile'), '--configfile', str(config_path),
                      '--cores', '2', '--quiet', 'all', '--', *targets])
        for name in ['First', 'Second']:
            (w / f'vtest/{name}/fastq_map.txt').write_text('reads.fq\trat_A\n')
        (w / 'reads').mkdir()
        (w / 'reads/reads.fq').write_text('')
        (w / 'ref/annotation.gtf').write_text('')
        (w / 'ref/star_index_100').mkdir()
        (w / 'ref/star_index_100/SAindex').write_text('')
        self.command(['snakemake', '-s', str(ROOT / 'Snakefile'), '--configfile', str(config_path),
                      '--profile', 'slurm', '--jobs', '2', '--dry-run', '--',
                      'vtest/First/star_out/rat_A.Aligned.sortedByCoord.out.bam',
                      'vtest/First/geno.bed', 'geno/individual/test/rat_A.vcf.gz',
                      'vtest/First/qc/test_snps.vcf.gz'])
        shared = pd.read_csv(w / 'geno/test.psam', sep=r'\s+')
        self.assertEqual(list(shared['#IID']), ['rat_A', 'rat_B', 'rat_C'])
        with pysam.VariantFile(str(w / 'geno/test.vcf.gz')) as vcf:
            self.assertEqual(list(vcf.header.samples), ['rat_A', 'rat_B', 'rat_C'])
            rows = {r.id: r for r in vcf}
        self.assertEqual(set(rows), {f'chr{c}:{p}' for c in [1, 2] for p in [10, 20, 60]})
        self.assertEqual(rows['chr1:20'].ref, 'A')
        self.assertNotIn('EAF', rows['chr1:20'].info)
        self.assertEqual(rows['chr1:20'].samples['rat_A']['GT'], (1, 1))
        bim = pd.read_csv(w / 'vtest/First/geno.bim', sep=r'\s+', header=None)
        self.assertEqual(set(bim[1]), {'chr1:10', 'chr1:20', 'chr2:10', 'chr2:20'})
        second_bim = pd.read_csv(w / 'vtest/Second/geno.bim', sep=r'\s+', header=None)
        self.assertEqual(set(second_bim[1]), set(rows))
        for name, expected in [('First', ['rat_A', 'rat_B']), ('Second', ['rat_A', 'rat_B', 'rat_C'])]:
            fam = pd.read_csv(w / f'vtest/{name}/geno.fam', sep=r'\s+', header=None)
            self.assertEqual(list(fam[1]), expected)
        with gzip.open(w / 'geno/test.alleles.tsv.gz', 'rt') as handle:
            alleles = dict((variant, (ref, alt)) for variant, ref, alt in (line.strip().split('\t') for line in handle))
        self.assertEqual(alleles['chr1:20'], ('A', 'G'))
        counts = json.loads((w / 'geno/test.work/chr1/variants.summary.json').read_text())['counts']
        self.assertEqual(counts, dict(keep=2, swap_ref_alt=1, exclude_reference_mismatch=1,
                                     exclude_duplicate_position=2, exclude_non_biallelic_snp=1))
        with pysam.VariantFile(str(w / 'geno/test.all_rats_exons.vcf.gz')) as vcf:
            self.assertEqual(list(vcf.header.samples), ids)
            self.assertEqual([r.id for r in vcf], ['chr1:10', 'chr1:20', 'chr2:10', 'chr2:20'])
        # PGEN retains dosages even though the BED and consumer VCF paths use hardcalls.
        self.command(['plink2', '--pfile', 'geno/test', '--export', 'vcf', 'vcf-dosage=DS-force', 'id-paste=iid', '--out', 'dosages'])
        with pysam.VariantFile(str(w / 'dosages.vcf')) as vcf:
            first = next(vcf)
            self.assertAlmostEqual(first.samples['rat_C']['DS'][0], 0.8, places=3)

    def command(self, command):
        env = dict(os.environ, XDG_CACHE_HOME=str(self.work / 'cache'))
        result = subprocess.run(command, cwd=self.work, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if result.returncode:
            self.fail(f'Command failed: {command}\n{result.stdout}')


if __name__ == '__main__':
    unittest.main()
