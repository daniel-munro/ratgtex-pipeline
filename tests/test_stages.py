"""Scientific workflow interfaces: tissue identity, project isolation, stage dispatch."""
import sys
import os
import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_stage
from workflow_config import datasets, pheast_config

ROOT = Path(__file__).resolve().parents[1]


class StageTests(unittest.TestCase):
    def setUp(self):
        (ROOT / '.snakemake').mkdir(exist_ok=True)
        self.temp = tempfile.TemporaryDirectory(dir=ROOT / '.snakemake')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.config = yaml.safe_load((ROOT / 'config/config.yaml').read_text())
        self.config.update(version='v5', tissues={'Combined': ['First', 'Second']},
                           datasets={'First': {'read_length': 100, 'fastq_path': 'reads/one'},
                                    'Second': {'read_length': 150, 'fastq_path': 'reads/two'}})
        source = self.root / 'Pantry'
        self.config['pantry']['source'] = str(source)
        for module in ['phenotyping', 'pheast']:
            (source / module).mkdir(parents=True)
            (source / module / 'Snakefile').write_text("configfile: 'config.yml'\n")
        for name, animal in [('First', 'rat1'), ('Second', 'rat2')]:
            project = self.root / 'v5' / name
            project.mkdir(parents=True)
            (project / 'rat_ids.txt').write_text(animal + '\n')
            (project / 'fastq_map.txt').write_text(f'{animal}.fq.gz\t{animal}\n')

    def setup_projects(self):
        with patch.object(run_stage.subprocess, 'check_output', return_value='revision\n'):
            run_stage.setup(self.config, self.root)
            run_stage.pantry_setup(self.config, self.root)

    def test_setup_merges_maps_and_reuses_sample_products(self):
        self.setup_projects()
        project = self.root / 'v5/Combined'
        self.assertEqual((project / 'rat_ids.txt').read_text(), 'rat1\nrat2\n')
        self.assertIn(str(self.root / 'reads/two/rat2.fq.gz'), (project / 'fastq_map.txt').read_text())
        self.assertFalse((project / 'phenos/intermediate').exists())
        run_stage.pantry_merge_setup(self.config, self.root)
        link = project / 'phenos/intermediate/expression/rat2'
        self.assertEqual(link.readlink(), self.root / 'v5/Second/phenos/intermediate/expression/rat2')
        for suffix in ['', '.bai']:
            bam = project / f'phenos/intermediate/bam/rat2.bam{suffix}'
            self.assertEqual(bam.readlink(), self.root / f'v5/Second/phenos/intermediate/bam/rat2.bam{suffix}')
        config = yaml.safe_load((project / 'phenos/config.yml').read_text())
        self.assertNotIn('paired_end', config)
        self.assertEqual(config['fastq_dir'], '/')
        self.assertEqual(datasets(self.config), ['First', 'Second'])
        qtl = yaml.safe_load((project / 'pheast/config.qtl.yml').read_text())
        twas = yaml.safe_load((project / 'pheast/config.twas.yml').read_text())
        self.assertIn('cross_modality', qtl['modalities'])
        self.assertNotIn('cross_modality', twas['modalities'])
        self.assertEqual(list(qtl['analyses']), ['qtl'])
        self.assertEqual(list(twas['analyses']), ['twas'])
        self.assertIn("configfile: 'config.twas.yml'", (project / 'pheast/Snakefile.twas').read_text())
        with self.assertRaises(FileExistsError):
            self.setup_projects()

    def test_merge_setup_preview_refresh_and_conflicts(self):
        self.setup_projects()
        args = SimpleNamespace(dry_run=True, plan=False)
        with patch.object(run_stage, 'ROOT', self.root):
            run_stage.run(self.config, 'pantry_merge_setup', self.root / 'config/config.yaml', args)
        destination = self.root / 'v5/Combined/phenos/intermediate/expression/rat2'
        self.assertFalse(destination.is_symlink())
        # Source products become available only after source phenotyping.
        source = self.root / 'v5/Second/phenos/intermediate/expression/rat2'
        source.mkdir(parents=True)
        os.utime(source, (100, 100))
        run_stage.pantry_merge_setup(self.config, self.root)
        os.utime(destination, (200, 200), follow_symlinks=False)
        run_stage.pantry_merge_setup(self.config, self.root)
        self.assertGreater(destination.lstat().st_mtime, 200)
        self.assertEqual(source.stat().st_mtime, 100)
        destination.unlink()
        destination.symlink_to(self.root / 'wrong')
        with self.assertRaisesRegex(ValueError, 'existing link'):
            run_stage.pantry_merge_setup(self.config, self.root)

    def test_combined_modality_headers_and_feature_groups(self):
        output = self.root / 'phenos/output'
        output.mkdir(parents=True)
        for modality, feature in [('expression', 'geneA'), ('splicing', 'geneB__intron')]:
            with gzip.open(output / f'{modality}.bed.gz', 'wt') as handle:
                handle.write(f'#chr\tstart\tend\tphenotype_id\trat1\nchr1\t20\t21\t{feature}\t0.5\n')
        command = [sys.executable, str(ROOT / 'scripts/combine_modalities.py'), str(self.root), 'expression', 'splicing']
        subprocess.run(command, check=True)
        self.assertEqual((output / 'cross_modality.phenotype_groups.txt').read_text(),
                         'expression:geneA\tgeneA\nsplicing:geneB__intron\tgeneB\n')
        with gzip.open(output / 'cross_modality.bed.gz', 'rt') as handle:
            self.assertEqual(len(handle.readlines()), 3)
        with gzip.open(output / 'splicing.bed.gz', 'wt') as handle:
            handle.write('#chr\tstart\tend\tphenotype_id\tDIFFERENT_ANIMAL\n')
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Sample header differs', result.stderr)

    def test_duplicate_animals_require_a_decision(self):
        (self.root / 'v5/Second/rat_ids.txt').write_text('rat1\n')
        with self.assertRaisesRegex(ValueError, 'overlapping animals'):
            self.setup_projects()
        self.assertFalse((self.root / 'v5/First/phenos').exists())

    def test_processing_stages_only_dispatch_to_ratgtex(self):
        args = SimpleNamespace(profile='slurm', jobs=10, dry_run=True, plan=False)
        for stage in ['genotypes', 'align', 'identity_qc', 'cross_modality', 'qtl_extra']:
            with self.subTest(stage=stage), patch.object(run_stage, 'ROOT', self.root), patch.object(run_stage.subprocess, 'run') as mock:
                run_stage.run(self.config, stage, self.root / 'config/config.yaml', args)
                mock.assert_called_once()
                self.assertEqual(mock.call_args.kwargs['cwd'], self.root)
                command = mock.call_args.args[0]
                self.assertEqual(command[1], stage)
                self.assertIn('--dry-run', command)
                self.assertIn('slurm', command)
                self.assertEqual(command[command.index('--configfile') + 1], str(self.root / 'config/config.yaml'))

    def test_general_setup_exports_lists_without_pantry_or_processing(self):
        directory = self.root / 'config'
        directory.mkdir()
        self.config['tissues'] = {'Combined': ['First', 'Second'], 'First': ['First']}
        self.config['datasets']['Unused'] = dict(self.config['datasets']['First'])
        args = SimpleNamespace(profile='slurm', jobs=10, dry_run=True, plan=False)
        with patch.object(run_stage, 'ROOT', self.root), patch.object(run_stage.subprocess, 'run') as mock:
            run_stage.run(self.config, 'setup', directory / 'config.yaml', args)
            self.assertFalse((directory / 'datasets.txt').exists())
            self.assertFalse((self.root / 'v5/Combined').exists())
            args.dry_run = False
            # General setup does not need a Pantry checkout or git metadata.
            with patch.object(run_stage.subprocess, 'check_output', side_effect=AssertionError('Pantry read')):
                run_stage.run(self.config, 'setup', directory / 'config.yaml', args)
            mock.assert_not_called()
        self.assertEqual((directory / 'datasets.txt').read_text(), 'First\nSecond\n')
        self.assertEqual((directory / 'tissues.txt').read_text(), 'Combined\nFirst\n')
        self.assertEqual((self.root / 'v5/Combined/rat_ids.txt').read_text(), 'rat1\nrat2\n')
        self.assertFalse((self.root / 'v5/First/phenos').exists())
        # A later identity correction refreshes combined inputs without touching Pantry.
        (self.root / 'v5/First/rat_ids.txt').write_text('corrected\n')
        (self.root / 'v5/First/fastq_map.txt').write_text('rat1.fq.gz\tcorrected\n')
        (self.root / 'v5/First/phenos').mkdir()
        sentinel = self.root / 'v5/First/phenos/config.yml'
        sentinel.write_text('keep this config\n')
        run_stage.setup(self.config, self.root)
        self.assertEqual((self.root / 'v5/Combined/rat_ids.txt').read_text(), 'corrected\nrat2\n')
        self.assertIn('corrected', (self.root / 'v5/Combined/fastq_map.txt').read_text())
        self.assertEqual(sentinel.read_text(), 'keep this config\n')

    def test_pantry_setup_is_separate_and_does_not_launch_processing(self):
        run_stage.setup(self.config, self.root)
        args = SimpleNamespace(profile='slurm', jobs=10, dry_run=True, plan=False)
        path = self.root / 'config/config.yaml'
        with patch.object(run_stage, 'ROOT', self.root), patch.object(run_stage.subprocess, 'run') as mock:
            run_stage.run(self.config, 'pantry_setup', path, args)
            self.assertFalse((self.root / 'v5/First/phenos').exists())
            args.dry_run = False
            with patch.object(run_stage.subprocess, 'check_output', return_value='revision\n'):
                run_stage.run(self.config, 'pantry_setup', path, args)
            mock.assert_not_called()
        settings = yaml.safe_load((self.root / 'v5/First/phenos/config.yml').read_text())
        self.assertEqual(settings['ref_genome'], str(self.root / self.config['ref_genome']))

    def test_pantry_setup_rejects_stale_combined_membership(self):
        run_stage.setup(self.config, self.root)
        (self.root / 'v5/Second/rat_ids.txt').write_text('corrected\n')
        with patch.object(run_stage.subprocess, 'check_output', return_value='revision\n'):
            with self.assertRaisesRegex(ValueError, 'rerun setup'):
                run_stage.pantry_setup(self.config, self.root)
        self.assertFalse((self.root / 'v5/First/phenos').exists())

    def test_pantry_execution_is_not_a_runner_stage(self):
        args = SimpleNamespace(profile='slurm', jobs=10, dry_run=False, plan=False)
        for stage in ['phenotypes', 'qtl', 'twas']:
            with self.subTest(stage=stage), patch.object(run_stage.subprocess, 'run') as mock:
                with self.assertRaisesRegex(ValueError, 'Run Pantry.*directly'):
                    run_stage.run(self.config, stage, self.root / 'config/config.yaml', args)
                mock.assert_not_called()


if __name__ == '__main__':
    unittest.main()
