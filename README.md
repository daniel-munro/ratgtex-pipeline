# RatGTEx Pipeline

This is the code used to process data for the [RatGTEx Portal](https://ratgtex.org). It is built on [Snakemake](https://snakemake.github.io/), a Python-based framework for reproducible data analysis. Since things might not work the first time with a new dataset, it makes it easy to process data iteratively. You can run part of the pipeline on a subset of the data, and then run it on the full dataset without regenerating existing files.

The main steps of the pipeline are:

1. Align RNA-Seq reads using [STAR](https://github.com/alexdobin/STAR).
2. Checking for and fixing sample mixups.
3. Quantify RNA phenotypes using [Pantry](https://github.com/PejLab/Pantry).
4. Map cis-eQTLs and trans-eQTLs using [tensorQTL](https://github.com/broadinstitute/tensorqtl) in various modes. Some of these are built into the Pantry Pheast module, while others are run in the main RatGTEx snakemake pipeline.
5. Compute TWAS weights for all RNA phenotypes using the Pantry Pheast module. These are used for the Rat TWAS Hub.
6. Calculate cis-eQTL effect size (allelic fold change) using [aFC.py](https://github.com/secastel/aFC).

## Setup

### Software environment

RatGTEx orchestration requires Python with PyYAML and Snakemake with the Slurm executor plugin. Alignment/QC and supplementary QTL rules also use the tools invoked in `steps/` (including STAR, samtools, bcftools, PLINK, GATK, tensorQTL, pandas, pysam, NumPy, fastparquet, bgzip, and tabix). Use your installed scientific software environment for those rules. Pantry supplies its own environment specifications; see its documentation for phenotyping and Pheast dependencies. There is no root RatGTEx `environment.yml`.

### Other software

To get `aFC.py`:

```
cd tools
git clone git@github.com:secastel/aFC.git
```

### Snakemake profile

When you run snakemake, you specify a profile that determines how steps get run. Here is the config file I use on a computing cluster with slurm scheduling:

`~/.config/snakemake/slurm/config.yaml`:

```yaml
executor: slurm
default-resources:
  runtime: "4h"
  mem_mb: 8000
  slurm_account: "our-account"
  slurm_partition: "our-cpu-partition"
  slurm_extra: "'--mail-type=FAIL --mail-user=my@email.address'"
set-resources:
  # `tasks_per_gpu: 0` prevents the error `srun: fatal: SLURM_TRES_PER_TASK is mutually exclusive with --ntasks-per-gpu and SLURM_NTASKS_PER_GPU`
  tensorqtl_cis:
    slurm_partition: "our-gpu-partition"
    gres: "gpu:1"
    tasks_per_gpu: 0
  tensorqtl_cis_independent:
    slurm_partition: "our-gpu-partition"
    gres: "gpu:1"
    tasks_per_gpu: 0
  tensorqtl_trans:
    slurm_partition: "our-gpu-partition"
    gres: "gpu:1"
    tasks_per_gpu: 0
  tensorqtl_cis_nominal:
    slurm_partition: "our-gpu-partition"
    gres: "gpu:1"
    tasks_per_gpu: 0
use-conda: true
latency-wait: 60
```

This uses the `snakemake-executor-plugin-slurm` Snakemake plugin. The `tensorqtl` steps should be run on GPU for reasonable runtime, so they are specified by name here to override the default resources. Adjust as needed for your cluster. Additional resources are specified within some of the snakemake rules, which are passed to slurm when those jobs are submitted. Alternatively, you can run snakemake on an interactive node.

#### `config/config.yaml`

`config/config.yaml` defines source **datasets**, analysis **tissues**, the shared genotype panel, modalities, reference files, and Pantry settings.

- A **dataset** is a tissue-specific sequencing collection from a source study, such as NAcc1.
- A **tissue** is the tissue-level analysis assembled from one or more datasets, such as NAcc from NAcc1–4 or Brain from Brain alone.
- The **genotype panel** is the shared round 11.3 PGEN fileset selected by `genotypes.panel`. Its indexed VCF serves alignment, QC, and aFC; tissue PLINK files serve Pantry.
- A **study** describes the originating project, potentially covering multiple tissues (e.g. IL_LHb_NAcc_OFC_PL). A **cohort** describes a defined group of animals within a study. Neither is currently a separate configuration entity.

The `tissues` mapping explicitly connects each tissue to its constituent datasets.

```yaml
version: "v5"
default_stage: align
tissues:
  IL: [IL]
  NAcc: [NAcc1, NAcc2, NAcc3, NAcc4]
datasets:
  IL:
    read_length: 100
    fastq_path: "fastq/IL_LHb_NAcc_OFC_PL"
    paired_end: false
  # Define each NAcc source dataset here too.
```

Only datasets referenced by `tissues` are selected. Tissue membership replaces the old `run`/`merged_tissues` lists and numeric-suffix convention. `modalities` defines each phenotype's Pantry producer and grouping; `cross_modality: true` combines all configured modalities for grouped cis-QTL analysis. TWAS uses the individual modalities because current Pheast restricts cross-modality analysis to QTL mapping.

`pantry.source` identifies the Pantry checkout copied during `pantry_setup`. `pantry.phenotyping` and `pantry.pheast` supply settings for generated project configs; `pantry.merge_products` declares the reusable sample-level products for each producer. Optional producer settings can be supplied under `pantry.phenotyping.modality_groups.<producer>`; output file lists are generated from `modalities`. Parameters still inside Pantry rules are not moved by this refactor.

The checked-in configuration selects v5 and round 11.3 genotypes for all datasets. Prepare the v5 dataset sample lists and FASTQ maps before setup. Paths in the main config are relative to the RatGTEx repository; generated Pantry configs use absolute input paths.

The current release configuration lives in `config/`:

```text
config/
    config.yaml       # Authoritative processing configuration
    tissue_info.tsv   # Tissue names, anatomical categories, and colors
    datasets.txt      # Generated list of selected source datasets
    tissues.txt       # Generated list of tissue-level analyses
```

Paths inside `config/config.yaml` remain relative to the repository root. Tissue metadata is maintained separately; `setup` generates the two ID lists from the main config while assembling combined-tissue sample lists/maps:

```shell
python scripts/run_stage.py setup
```

With `--config <path>`, lists are written beside that configuration file. `scripts/qc/sample_tracking.py` defaults to `config/datasets.txt` (its existing `--tissues` argument can select a different dataset-list file). Archived v4 configuration and lists remain under `v4/`; current code uses `config/`.

### Dataset-specific input files

#### FASTQ files

The RNA-seq FASTQ files can be single-read or paired-end. There can be more than one file or file pair per sample, and all will be used.

#### `{version}/{dataset}/fastq_map.txt`

A tab-delimited file with no header containing the paths to each FASTQ file and the rat IDs they correspond to. Or, for paired-end reads, each row contains the first FASTQ path, second FASTQ path, and rat ID per file pair.

- If multiple files map to the same ID, i.e. the ID appears in multiple rows, reads from those files will be aligned into one BAM file.
- You can use the `fastq_path` parameter in `config/config.yaml` to specify the encompassing directory as an absolute or relative path. That way `fastq_map.txt` can just contain the remainder of the path to each file (including any subdirectories as necessary).
- Any listed files whose rat IDs are not in `{version}/{dataset}/rat_ids.txt` will be ignored.
- If a dataset includes a mix of single-end and paired-end samples, omit `paired_end` from the dataset config so the pipeline infers pairing per row.

#### `{version}/{dataset}/rat_ids.txt`

A file listing the rat IDs for the dataset, one per line. This list determines which samples are included in the processing.

### v5 genotypes

Round 11.3 inputs are chromosome-specific PGEN/PVAR/PSAM files. `genotypes` in the main config supplies their prefix template, the output panel name, and an explicit chromosome list (autosomes 1–20). All datasets use this shared panel; there are no dataset-specific genotype panels in v5.

After general `setup` has assembled the tissue sample lists, run the `genotypes` stage. It uses the main Snakefile and reads sample lists without loading FASTQ maps or Pantry projects:

```shell
python scripts/run_stage.py genotypes -n
# Remove -n to prepare the shared panel, consumer VCF, and tissue BED files.
```

The workflow maps RFIDs to round 11.3 sample IDs, subsets each chromosome to the union of requested rats, corrects REF/ALT orientation against the configured FASTA, and merges the chromosomes. RFIDs come from the CSV configured by `genotypes.metadata`. Metadata `library_name`, `rfid`, and `barcode` fields construct the two PSAM join-key forms (`library_name_rfid` and `library_name_rfid_barcode`); genotype IDs are not parsed. Every PSAM sample must map to one metadata RFID, and conflicting assignments cause failure. Metadata-only rats absent from the filtered PSAM are not treated as available genotypes. Missing or ambiguous requested rats cause failure; no replicate is selected implicitly. Every chromosome's PSAM must match, including sample order and sex.

The shared outputs are `geno/ratgtex_v5_round11_3.{pgen,pvar,psam}`, an indexed `.vcf.gz` for STAR alignment, identity QC, and aFC, and `.alleles.tsv.gz` with variant ID, REF, and ALT (no header). PGEN retains available dosage and phase; consumer VCFs export hardcalls. The `.work/` directory holds the original-ID/RFID keep and rename tables, per-dataset sample coverage/counts, intermediate chromosome panels, and per-chromosome variant audit/summary files. Coverage reports are Snakemake logs so they remain available when identity validation fails.

Only biallelic A/C/G/T SNPs that match the FASTA in either orientation are retained. SNPs matching only ALT are flipped. Unresolved reference mismatches and all records at duplicate positions are excluded and audited. Excluding every duplicate-position record avoids choosing an arbitrary allele at a site. Missing or repeated input variant IDs cause failure. Final IDs are `chrN:position`, with chromosome labels `chr1`–`chr20`. Upstream INFO fields are dropped because cohort frequencies and allele-specific annotations would be stale after sample subsetting or allele swaps. No allele-frequency filter is applied to the shared panel.

Tissue PLINK files are derived directly from the shared PGEN panel, using the tissue's sample list and `--maf 0.01 --mac 2 --max-alleles 2`. Filtering runs on the extracted hardcalls, preserving the previous GT-based frequency/count semantics even when PGEN contains dosages. The `genotypes` stage requests these files for every configured analysis tissue. Individual tissues can also be requested directly:

```shell
snakemake v5/NAcc/geno.bed v5/PL/geno.bed --profile slurm -j10 -n
```

For unexpected-identity searches, prepare an exon-only VCF containing **all original round 11.3 samples**, with the same reference validation and variant IDs:

```shell
snakemake geno/ratgtex_v5_round11_3.all_rats_exons.vcf.gz \
    geno/ratgtex_v5_round11_3.all_rats_exons.snps.txt --profile slurm -j10 -n
```

This optional target extracts exons before exporting each chromosome, then concatenates chromosomes in configured order. It retains original genotyping IDs so unexpected rats remain identifiable. The QC comparison rule uses this release-specific VCF and SNP list. Large PLINK jobs use a project-local `tmpdir` resource. Archived v1–v4 preparation scripts remain historical documentation.

## Running

Start by supplying `{version}/{dataset}/rat_ids.txt` and `fastq_map.txt` for every selected **source dataset**. For existing datasets, copy these from `samples/rat_ids_current/rat_ids-{dataset}.txt` and `samples/fastq_map_current/fastq_map-{dataset}.txt`. The current list/map files are the source of the derived sample tracking tables, not outputs generated from those tables. New datasets require list/map inputs derived from their sequencing metadata; the v5 VMH inputs are documented in [v5/VMH/README.md](v5/VMH/README.md). For a tissue containing a single source dataset, those files also serve as the tissue's files.

Then proceed in the order below. General setup assembles **combined-tissue** lists/maps (e.g. NAcc from NAcc1–4); alignment uses the source-dataset inputs. Pantry projects are created after identity review.

| Step | Stage or manual action | Result / prerequisite |
| ---- | ---------------------- | --------------------- |
| 1 | `setup` | Assemble combined-tissue sample lists/maps and export `config/datasets.txt` and `tissues.txt`. Reads the manually supplied source files. |
| 2 | `genotypes` | Prepare the shared round 11.3 PGEN panel, indexed consumer VCF, allele table, and BED/BIM/FAM files for all analysis tissues. |
| 3 | `align` | Align source datasets using their FASTQ maps and individual genotype VCFs; produce BAM indexes and STAR summaries. |
| 4 | `identity_qc` | Compare RNA and genotype identities for source datasets. |
| 5 | Manual identity review | Resolve relabeling/removal decisions. If inputs change, rerun `setup` and `genotypes`, rebuild affected alignments, and repeat QC. |
| 6 | `pantry_setup` | Copy Pantry projects, write initial configs, and link prepared inputs/sample products. Uses the reviewed sample membership. |
| 7 | Manual Pantry phenotyping | Run source-dataset projects, then combined-tissue projects after their constituent datasets finish. |
| 8 | `cross_modality` | Combine tissue phenotype tables after Pantry phenotyping. |
| 9 | Manual Pheast QTL / TWAS | Run each tissue's Pheast project; QTL uses cross-modality phenotypes when enabled. TWAS uses individual modalities. |
| 10 | `qtl_extra` | Generate significant-pair, expression p-value, aFC, and trans-QTL outputs after Pheast QTL results/covariates are ready. Requires expression. |

From the RatGTEx repository:

```shell
python scripts/run_stage.py setup
python scripts/run_stage.py genotypes -n
python scripts/run_stage.py align -n
python scripts/run_stage.py identity_qc -n
# Review identities and repeat affected earlier steps before proceeding.
python scripts/run_stage.py pantry_setup -n
# Run Pantry phenotyping manually, as described below.
python scripts/run_stage.py cross_modality -n
# Run Pheast manually, as described below.
python scripts/run_stage.py qtl_extra -n
```

`--profile slurm` and `--jobs 10` are defaults for processing stages. Remove `-n` to execute a stage. `setup -n` and `pantry_setup -n` preview their Python actions without writing files; processing stages pass `-n` to Snakemake and require their input files. `--plan` previews actions/commands without executing them. The runner does not automatically run earlier stages or launch Pantry analyses. Snakemake resolves missing upstream processing products within the requested stage's DAG.

General `setup` performs no genotype processing or Pantry copying. It can be rerun after updating source membership/maps to refresh combined-tissue inputs and the configuration ID lists. It leaves source files and existing Pantry projects alone. `pantry_setup` performs no processing and refuses existing Pantry project directories, preserving their configs/results.

All processing stages use the main RatGTEx Snakefile; their targets are available directly:

```shell
snakemake genotypes --profile slurm -j10 -n
snakemake align --profile slurm -j10 -n
snakemake identity_qc --profile slurm -j10 -n
snakemake cross_modality --profile slurm -j10 -n
snakemake qtl_extra --profile slurm -j10 -n
```

Without a target, Snakemake uses `default_stage` (initially `align`). Its supported values are `genotypes`, `align`, `identity_qc`, and `qtl_extra`. To process a subset, use a config containing only the desired tissues, or request individual output paths directly. FASTQ maps are loaded only by alignment rules, so a genotype target does not require them.

### Identity review before Pantry setup

Review identity matches manually. For unresolved samples, list their IDs plus a positive control in `{version}/{dataset}/qc/samples_without_matches.txt` and request `{version}/{dataset}/qc/all_rats_summary.tsv` with Snakemake. This builds the optional all-rats exon VCF as needed.

Record relabeling/removal decisions, update the current rat-ID lists and FASTQ maps and their release-local copies, then regenerate the derived sample tracking tables. Run `setup` to refresh combined-tissue files, then `genotypes` to refresh genotype inputs. Delete BAMs for changed identities and regenerate alignments and dependent outputs. Genotype-aware alignment makes rebuilding these BAMs necessary. Repeat identity QC before `pantry_setup`.

After Pantry setup, source-dataset phenotyping projects can run concurrently. A combined tissue must wait for its constituent source datasets' sample products; unrelated tissues can proceed independently. Run RatGTEx `cross_modality` after tissue phenotyping, Pheast QTL/TWAS manually, then RatGTEx `qtl_extra` after the relevant Pheast QTL results and covariates are ready.

The aggregate targets cover all configured tissues. To proceed with one ready tissue while others are still running, request its output files directly:

```shell
snakemake v5/Brain/phenos/output/cross_modality.bed.gz --profile slurm -j10 -n
snakemake v5/Brain/geno.bed --profile slurm -j10 -n
snakemake v5/Brain/Brain.expression.cis_qtl_signif.txt.gz --profile slurm -j10 -n
```

These RatGTEx targets consume existing Pantry products and do not launch Pantry to generate missing products. To force a tissue genotype rebuild, add `--forcerun pgen_to_plink` to its output-target command; normally Snakemake tracks sample-list changes automatically.

### Pantry projects and tissue assembly

`pantry_setup` creates `{version}/{name}/phenos` and `pheast` projects for selected source datasets and analysis tissues. It refuses existing Pantry project directories, protecting previous code and results; do not use it to migrate the existing v4 copies. It records the Pantry revision, working-tree status, and main configuration in `ratgtex_setup.json`.

Run phenotyping directly inside the project you want to process:

```shell
cd v5/Brain/phenos
snakemake --profile slurm -j20 -n
```

Remove `-n` to run. Start another tissue in a separate session to run it concurrently. For merged NAcc, run the same command inside `v5/NAcc/phenos` after NAcc1–4 have finished.

Each Pheast project has `config.qtl.yml`/`Snakefile.qtl` and `config.twas.yml`/`Snakefile.twas`. Select the desired entry point yourself. Initial `config.yml` defaults to QTL and can be edited for additional interactive analyses. For example:

```shell
cd v5/Brain/pheast
snakemake --snakefile Snakefile.qtl --profile slurm -j10 -n
# Or, for TWAS:
snakemake --snakefile Snakefile.twas --profile slurm -j10 -n
```

`pantry_setup` writes the initial Pantry configs once. If parameters change later, edit the relevant project configs manually; changing the main config does not rewrite existing Pantry projects. If sample membership changes after Pantry setup, rerun general `setup` and `genotypes`, and update affected Pantry links/configs manually; `pantry_setup` does not overwrite existing projects.

Combined tissues reuse source-dataset sample-level products through symlinks, then compute tissue-dependent products in their own projects. Their maps contain absolute FASTQ paths, so mixed read lengths/layouts retain their original source quantification. Overlapping animal IDs across constituent datasets raise an error rather than silently deduplicating libraries. Decide how those libraries should be combined before proceeding. The old version-specific setup scripts remain historical documentation; the stage runner replaces them for new project setup.

### Workflow checks

```shell
python -m unittest discover -s tests -v
```

These checks exercise setup isolation, sample identity handling, config generation, and stage dispatch. Use Snakemake dry runs to validate the actual selected data and dependencies.

## Help

Email me or file an issue on this GitHub repository if anything isn't working, and I'll be happy to fix it.
