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
from genotype_samples import metadata_sample_ids, prepare_samples


class GenotypeTests(unittest.TestCase):
    def setUp(self):
        (ROOT / '.snakemake').mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / '.snakemake')
        self.addCleanup(self.temp.cleanup)
        self.work = Path(self.temp.name)
        self.metadata = self.work / 'metadata.csv'
        pd.DataFrame([
            ('rat_A', 'Rattaca1', 'x'), ('rat_A', 'Riptide2', 'x'),
            ('00079299E6', 'ddGBS_HSR_96_WF1', 'AATATGG'),
            ('absent', 'arbitrary_library', 'barcode'),
        ], columns=['rfid', 'library_name', 'barcode']).to_csv(self.metadata, index=False)

    def test_sample_coverage_and_ambiguity(self):
        w = self.work
        sample = w / 'rats.txt'
        sample.write_text('rat_A\nabsent\n')
        psam = w / 'samples.psam'
        psam.write_text('#IID\tSEX\nRattaca1_rat_A_x\t1\n')
        outputs = [w / name for name in ['keep', 'rename', 'coverage', 'counts']]
        with self.assertRaisesRegex(ValueError, 'Missing or ambiguous'):
            prepare_samples([sample], ['dataset'], [psam], self.metadata, *outputs)
        coverage = pd.read_csv(outputs[2], sep='\t')
        self.assertEqual(list(coverage.status), ['present', 'missing'])
        sample.write_text('rat_A\n')
        psam.write_text(psam.read_text() + 'Riptide2_rat_A_x\t1\n')
        with self.assertRaisesRegex(ValueError, 'Missing or ambiguous'):
            prepare_samples([sample], ['dataset'], [psam], self.metadata, *outputs)
        self.assertEqual(pd.read_csv(outputs[2], sep='\t').status[0], 'ambiguous')
        other = w / 'other.psam'
        other.write_text('#IID\tSEX\nRattaca1_rat_A_x\t2\n')
        with self.assertRaisesRegex(ValueError, 'PSAM differs'):
            prepare_samples([sample], ['dataset'], [psam, other], self.metadata, *outputs)

    def test_ddgbs_wf1_library_component_is_not_part_of_rfid(self):
        w = self.work
        sample = w / 'rats.txt'
        sample.write_text('00079299E6\n')
        psam = w / 'samples.psam'
        psam.write_text('#IID\tSEX\nddGBS_HSR_96_WF1_00079299E6\t1\n')
        outputs = [w / name for name in ['keep', 'rename', 'coverage', 'counts']]
        prepare_samples([sample], ['Brain'], [psam], self.metadata, *outputs)
        coverage = pd.read_csv(outputs[2], sep='\t', dtype=str)
        self.assertEqual(coverage.rfid.tolist(), ['00079299E6'])
        self.assertEqual(coverage.status.tolist(), ['present'])
        self.assertEqual(outputs[1].read_text().splitlines()[-1],
                         '0\tddGBS_HSR_96_WF1_00079299E6\t0\t00079299E6')

    def test_metadata_join_supports_arbitrary_library_names_and_rejects_conflicts(self):
        w = self.work
        sample = w / 'rats.txt'
        sample.write_text('animal_with_underscores\n')
        psam = w / 'samples.psam'
        psam.write_text('#IID\tSEX\narbitrary_library_v3_animal_with_underscores_tag\t1\n')
        pd.DataFrame([('animal_with_underscores', 'arbitrary_library_v3', 'tag')],
                     columns=['rfid', 'library_name', 'barcode']).to_csv(self.metadata, index=False)
        outputs = [w / name for name in ['keep', 'rename', 'coverage', 'counts']]
        prepare_samples([sample], ['dataset'], [psam], self.metadata, *outputs)
        self.assertEqual(pd.read_csv(outputs[2], sep='\t').status.tolist(), ['present'])
        psam.write_text(psam.read_text() + 'unannotated_sample\t2\n')
        with self.assertRaisesRegex(ValueError, 'absent from genotype metadata'):
            prepare_samples([sample], ['dataset'], [psam], self.metadata, *outputs)
        pd.DataFrame([('B_C', 'A', 'x'), ('C', 'A_B', 'y')],
                     columns=['rfid', 'library_name', 'barcode']).to_csv(self.metadata, index=False)
        with self.assertRaisesRegex(ValueError, 'Conflicting metadata RFIDs'):
            metadata_sample_ids(self.metadata)

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
        pd.DataFrame([('rat_A', 'Rattaca1', 'x'), ('rat_B', 'ddGBS_HSR_1', 'barcode'),
                      ('rat_C', 'LSWrnaseq', 'LSWNNNNN'), ('unexpected', 'Riptide1', 'x')],
                     columns=['rfid', 'library_name', 'barcode']).to_csv(self.metadata, index=False)
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
        config = dict(version='vtest', genotypes=dict(panel='test', ld_reference_tissue='First', source_prefix='input/chr{chrom}', metadata=str(self.metadata), chromosomes=[1, 2]),
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
        ld_ids = []
        for chrom in [1, 2]:
            prefix = w / f'geno/test.LDREF/First_vtest.chr{chrom}'
            ld_bim = pd.read_csv(str(prefix) + '.bim', sep=r'\s+', header=None)
            self.assertEqual(set(ld_bim[0]), {chrom})
            self.assertEqual(list(ld_bim[1]), list(bim.loc[bim[0] == f'chr{chrom}', 1]))
            ld_fam = pd.read_csv(str(prefix) + '.fam', sep=r'\s+', header=None)
            self.assertEqual(list(ld_fam[1]), ['rat_A', 'rat_B'])
            ld_ids.extend(ld_bim[1])
        self.assertEqual((w / 'geno/test.LD_SNPs.txt').read_text().splitlines(), ld_ids)
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
