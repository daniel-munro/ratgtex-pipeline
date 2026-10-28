"""Prepare and download the two archival HS rat RNA-seq datasets from ENA."""

import argparse
import csv
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess


PROJECTS = {"PRJEB60407": ("Amyg", "SINGLE", 100),
            "PRJEB60349": ("Heart", "PAIRED", 101)}


def read_table(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def prepare(root):
    files, samples = [], []
    for project, (tissue, layout, length) in PROJECTS.items():
        runs = read_table(root / "metadata" / f"{project}.tsv")
        supplied = read_table(root / "metadata" / f"filereport_read_run_{project}.tsv")
        assert len({r['run_accession'] for r in runs}) == len(runs)
        assert {r['run_accession'] for r in runs} == {r['run_accession'] for r in supplied}
        original = {r['run_accession']: r for r in supplied}
        for row in runs:
            run = row['run_accession']
            assert row['study_accession'] == project
            assert row['library_strategy'] == 'RNA-Seq'
            assert row['library_layout'] == layout
            assert int(row['base_count']) == length * int(row['read_count'])
            assert row['fastq_ftp'] == original[run]['fastq_ftp']
            urls = row['fastq_ftp'].split(';')
            checksums = row['fastq_md5'].split(';')
            sizes = row['fastq_bytes'].split(';')
            assert len(urls) == len(checksums) == len(sizes) == (2 if layout == 'PAIRED' else 1)
            names = []
            for url, checksum, size in zip(urls, checksums, sizes):
                name = url.rsplit('/', 1)[1]
                names.append(name)
                assert len(checksum) == 32 and int(size) > 0
                files.append(dict(tissue=tissue, run_accession=run, filename=name,
                                  url='https://' + url, md5=checksum, bytes=size))
            samples.append(dict(tissue=tissue, study_accession=project,
                                run_accession=run, sample_accession=row['sample_accession'],
                                archive_sample=row['sample_alias'].split(':', 1)[1],
                                library_layout=layout, read_length=length,
                                fastq_1=names[0], fastq_2=names[1] if len(names) == 2 else ''))
    assert len({r['filename'] for r in files}) == len(files)
    for name, rows in [('download_manifest.tsv', files), ('sample_metadata.tsv', samples)]:
        with open(root / name, 'w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
    print(f"Prepared {len(samples)} runs, {len(files)} files, "
          f"{sum(int(r['bytes']) for r in files) / 1e12:.3f} TB compressed", flush=True)


def verify(path, row):
    if path.stat().st_size != int(row['bytes']):
        raise ValueError(f"Size mismatch: {path}")
    with open(path, 'rb') as handle:
        digest = hashlib.file_digest(handle, 'md5').hexdigest()
    if digest != row['md5']:
        raise ValueError(f"MD5 mismatch: {path}: expected {row['md5']}, observed {digest}")


def download(root, row):
    path = root / row['tissue'] / row['filename']
    # Download completion and checksum verification are separate stages.
    if path.exists():
        if path.stat().st_size != int(row['bytes']):
            raise ValueError(f"Size mismatch: {path}")
        print(f"Already downloaded (MD5 unchecked): {path.name}", flush=True)
        return
    else:
        partial = path.with_name(path.name + '.part')
        if not partial.exists() or partial.stat().st_size != int(row['bytes']):
            subprocess.run(['curl', '--fail', '--location', '--retry', '5',
                        '--retry-all-errors', '--connect-timeout', '60',
                        '--no-progress-meter', '--continue-at', '-',
                        '--output', str(partial), row['url']], check=True)
        if partial.stat().st_size != int(row['bytes']):
            raise ValueError(f"Size mismatch: {partial}")
        partial.rename(path)
    print(f"Downloaded (MD5 unchecked): {path.name}", flush=True)


def process(root, row, action):
    """Record each failure and continue; the overall command still fails."""
    try:
        if action == 'verify':
            verify(root / row['tissue'] / row['filename'], row)
            print(f"Verified {row['filename']}", flush=True)
        elif action == 'redownload':
            # Never resume suspect bytes or replace the original before verification.
            path = root / row['tissue'] / row['filename']
            fresh = path.with_name(path.name + '.fresh')
            subprocess.run(['curl', '--fail', '--location', '--connect-timeout', '60',
                            '--no-progress-meter', '--output', str(fresh), row['url']],
                           check=True)
            verify(fresh, row)
            fresh.replace(path)
            print(f"Replaced with verified fresh download: {path.name}", flush=True)
        else:
            download(root, row)
        return dict(filename=row['filename'], status='ok', error='')
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"FAILED {row['filename']}: {error}", flush=True)
        return dict(filename=row['filename'], status='failed', error=str(error))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'download', 'verify', 'redownload'])
    parser.add_argument('--directory', type=Path, default=Path('fastq/Amyg_Heart'))
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--failure-log', type=Path,
                        help='Verifier stdout log selecting files for redownload')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.workers < 1:
        parser.error('--workers must be positive')
    if args.action == 'redownload' and args.failure_log is None:
        parser.error('redownload requires --failure-log')
    if args.dry_run and args.action != 'redownload':
        parser.error('--dry-run is supported only for redownload')
    if args.action == 'prepare':
        prepare(args.directory)
    else:
        rows = read_table(args.directory / 'download_manifest.tsv')
        if args.action == 'redownload':
            with args.failure_log.open() as handle:
                failed_names = [line.split()[1].removesuffix(':')
                                for line in handle if line.startswith('FAILED ')]
            if not failed_names or len(failed_names) != len(set(failed_names)):
                raise ValueError('Expected a nonempty, unique list of FAILED files')
            unknown = set(failed_names) - {row['filename'] for row in rows}
            if unknown:
                raise ValueError(f'Failed files absent from manifest: {sorted(unknown)}')
            rows = [row for row in rows if row['filename'] in set(failed_names)]
            print(f"Selected {len(rows)} files, "
                  f"{sum(int(row['bytes']) for row in rows) / 1e9:.2f} GB", flush=True)
            if args.dry_run:
                for row in rows:
                    print(row['filename'])
                return
        if args.action in ('download', 'redownload'):
            for tissue, _, _ in PROJECTS.values():
                (args.directory / tissue).mkdir(exist_ok=True)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(lambda row: process(args.directory, row, args.action), rows))
        report = args.directory / f'{args.action}_results.tsv'
        with report.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=['filename', 'status', 'error'], delimiter='\t')
            writer.writeheader()
            writer.writerows(results)
        failed = sum(r['status'] == 'failed' for r in results)
        print(f"{args.action}: {len(results) - failed} succeeded, {failed} failed; {report}", flush=True)
        if failed:
            raise SystemExit(1)


if __name__ == '__main__':
    main()
