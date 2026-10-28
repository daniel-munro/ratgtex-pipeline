# Archival HS rat Amyg and Heart downloads

Run commands from the RatGTEx repository root. Prepared on 2026-10-20 from
the supplied `../rat-datasets/archival_hs_rats_baud/data/filereport_read_run_*.tsv`
and ENA run reports, with sample identities checked against ArrayExpress SDRFs.

| Dataset | ENA study | ArrayExpress | Runs | Layout | Read length | Compressed GB |
| --- | --- | --- | ---: | --- | ---: | ---: |
| Amyg | PRJEB60407 | E-MTAB-12693 | 205 | single | 100 | 252.22 |
| Heart | PRJEB60349 | E-MTAB-12701 | 192 | paired | 101 | 1242.17 |

The prepared manifests and source metadata are in `fastq/Amyg_Heart/`:

- `metadata/filereport_read_run_PRJEB*.tsv`: copies of the supplied reports.
- `metadata/PRJEB*.tsv`: ENA reports including aliases, layout, read/base counts,
  URLs, compressed sizes, and MD5 checksums.
- `metadata/E-MTAB-*.sdrf.txt`: ArrayExpress sample annotations.
- `sample_metadata.tsv`: one row per run, retaining biological archive aliases.
- `download_manifest.tsv`: 589 expected FASTQs with HTTPS URLs and checksums.

Download on a node with internet access and at least 1.50 TB free (decimal units):

```bash
python scripts/download/amyg_heart.py download --workers 4
```

For a Slurm submission (adjust partition/account/time to local requirements):

```bash
sbatch --job-name=ena_Amyg_Heart --cpus-per-task=4 --mem=2G --time=7-00:00:00 \
  --output=fastq/Amyg_Heart/download-%j.log \
  --wrap='python -u scripts/download/amyg_heart.py download --workers 4'
```

FASTQs are written beneath `fastq/Amyg_Heart/Amyg/` and
`fastq/Amyg_Heart/Heart/`, retaining ENA run filenames. The manifest tissue
column determines the subdirectory; metadata, manifests, and results remain
in `fastq/Amyg_Heart/`. Sample maps retain basenames, with each dataset’s
`fastq_path` selecting its tissue directory. Transfers resume `.part` files and retry all curl transfer errors
up to five times, including interrupted connections (exit code 56). Each
completed file must pass the compressed-size check before being renamed.
Download completion does not imply MD5 verification. Rerunning reuses files with
the expected size, including full-size partial files. Checksums are checked only
by the separate `verify` command below. Each pass attempts every manifest entry,
writes `download_results.tsv` or `verify_results.tsv`, and exits nonzero at the end
if any files failed. Verification reports expected and observed mismatching MD5s.
Inspect failed files before deciding which need fresh downloads. Run only one downloader
against this directory at a time. No SRA conversion or large temporary files are
needed; partial transfers stay beside their final outputs.

To verify every expected output independently:

```bash
python scripts/download/amyg_heart.py verify --workers 4
```

To replace failed files from a saved verifier stdout log, review the selection
and then submit from the repository root:

```bash
python scripts/download/amyg_heart.py redownload --failure-log amyg_heart_verify.txt --dry-run
sbatch scripts/download/amyg_heart_redownload.sbatch
```

The submission uses four workers, 8 GB memory, the `gdml` account and `cpu-gdml`
partition. Each selected file downloads to `<filename>.fresh` without resume or
automatic retries. Only after size and MD5 verification does it replace the
original. Failed fresh copies remain for inspection; rerunning starts those
downloads from scratch. Other selected files continue after a failure, and the
job exits nonzero if any fail. Results are saved to `redownload_results.tsv`,
preserving the earlier verification report. Do not run concurrent downloaders.
After this job succeeds, rerun the full verification command above.

To regenerate the manifests from the saved ENA reports:

```bash
python scripts/download/amyg_heart.py prepare
```

The ENA reports can be fetched again using this URL pattern for each study:

```text
https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB60407&result=read_run&fields=run_accession,study_accession,sample_accession,sample_alias,experiment_accession,library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes&format=tsv
```

SDRF source URLs:

- https://www.ebi.ac.uk/biostudies/files/E-MTAB-12693/E-MTAB-12693.sdrf.txt
- https://www.ebi.ac.uk/biostudies/files/E-MTAB-12701/E-MTAB-12701.sdrf.txt

## RatGTEx v5 integration

The curated `HS_rat_RNAseq_to_PLINK_ID_mapping_20261020.tsv` crosswalk was joined
one-to-one by ENA run accession, with public sample name, tissue, and study
checked independently. `plink_iid` is converted to `Outbred{plink_iid}`.
All 397 runs and rats within each tissue are unique. The following command
recreates the original unfiltered inputs and overwrites curated ID lists; do not
rerun it over the current exclusions or later identity-QC decisions:

```bash
python scripts/setup/amyg_heart.py ../rat-datasets/archival_hs_rats_baud/data/HS_rat_RNAseq_to_PLINK_ID_mapping_20261020.tsv
```

Authoritative lists/maps were written under `samples/rat_ids_current` and
`samples/fastq_map_current` and synchronized to `v5/Amyg` and `v5/Heart`.
Release-local `sample_metadata.tsv` files preserve the identity join and genotype
coverage. Both datasets are enabled in `config/config.yaml`, with tissue
annotations in `config/tissue_info.tsv`; tissue colors are unassigned.

All 20 source chromosome PSAMs agree. Coverage against the filtered round 11.3
panel is:

| Dataset | Included rats | Genotyped | Missing from filtered panel | Missing but in genotype metadata |
| --- | ---: | ---: | ---: | ---: |
| Amyg | 205 | 182 | 23 | 2 |
| Heart | 192 | 169 | 23 | 4 |

`fastq/Amyg_Heart/genotype_coverage.tsv` records each identity and distinguishes
filtered-panel availability from metadata presence. The 23 samples per tissue
missing from the filtered panel have been removed from both current and v5
`rat_ids.txt` lists, leaving 182 Amyg and 169 Heart samples. Exclusions are recorded
in `samples/removed.tsv`: 40 `not_genotyped` and six
`in_11_3_metadata_only` (two Amyg, four Heart). FASTQ maps, sample metadata,
and the download manifest retain every run for tracking and identity review.

After downloading, dry-run alignment with:

```bash
snakemake --profile slurm -n v5/Amyg/qc/Amyg.star_stats.tsv v5/Heart/qc/Heart.star_stats.tsv
```

Retain the manual identity-QC review before downstream analysis. All 589 FASTQs have been downloaded and verified.
