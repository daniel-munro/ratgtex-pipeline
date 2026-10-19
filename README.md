# RatGTEx Pipeline

This is the code used to process data for the [RatGTEx Portal](https://ratgtex.org). It is built on [Snakemake](https://snakemake.github.io/), a Python-based framework for reproducible data analysis. Since things might not work the first time with a new dataset, it makes it easy to process data iteratively. You can run part of the pipeline on a subset of the data, and then run it on the full dataset without regenerating existing files.

The main steps of the pipeline are:

1. Align RNA-Seq reads using [STAR](https://github.com/alexdobin/STAR).
2. Check for and fix sample mixups.
3. Quantify RNA phenotypes using [Pantry](https://github.com/PejLab/Pantry).
4. Map cis-eQTLs and trans-eQTLs using [tensorQTL](https://github.com/broadinstitute/tensorqtl) in various modes. Some of these are built into the Pantry Pheast module, while others are run in the main RatGTEx snakemake pipeline.
5. Compute TWAS weights for all RNA phenotypes using the Pantry Pheast module. These are used for the Rat TWAS Hub.
6. Calculate cis-eQTL effect size (allelic fold change) using [aFC.py](https://github.com/secastel/aFC).

## Setup

### Software environment

RatGTEx orchestration requires Python with PyYAML and Snakemake with the Slurm executor plugin. Alignment/QC and supplementary QTL rules also use the tools invoked in `steps/` (including STAR, samtools, bcftools, PLINK, GATK, tensorQTL, pandas, pysam, NumPy, fastparquet, bgzip, and tabix). Use your installed scientific software environment for those rules. Pantry supplies its own environment specifications; see its documentation for phenotyping and Pheast dependencies.

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

### Configuration

`config/config.yaml` defines source **datasets**, analysis **tissues**, the shared genotype panel, modalities, reference files, and Pantry settings.

- A **dataset** is a tissue-specific sequencing collection from a source study, such as NAcc1.
- A **tissue** is the tissue-level analysis assembled from one or more datasets, such as NAcc from NAcc1–4 or Brain from Brain alone.
- The **genotype panel** is the shared round 11.3 PGEN fileset selected by `genotypes.panel`. Its indexed VCF serves alignment, QC, and aFC; tissue PLINK files serve Pantry.

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

Only datasets referenced by `tissues` are selected. `modalities` defines each phenotype's Pantry producer and grouping; `cross_modality: true` combines all configured modalities for grouped cis-QTL analysis. TWAS uses the individual modalities because current Pheast restricts cross-modality analysis to QTL mapping.

`pantry.source` identifies the Pantry checkout copied during `pantry_setup`. `pantry.phenotyping` and `pantry.pheast` supply settings for generated project configs; `pantry.merge_products` declares the reusable sample-level products for each producer. Optional producer settings can be supplied under `pantry.phenotyping.modality_groups.<producer>`; output file lists are generated from `modalities`.

The checked-in configuration selects v5 and round 11.3 genotypes for all datasets. Prepare the v5 dataset sample lists and FASTQ maps before setup. Paths in the main config are relative to the RatGTEx repository; generated Pantry configs use absolute input paths.

The current release configuration lives in `config/`:

```text
config/
    config.yaml       # Authoritative processing configuration
    tissue_info.tsv   # Tissue names, anatomical categories, and colors
    datasets.txt      # Generated list of selected source datasets
    tissues.txt       # Generated list of tissue-level analyses
```

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

## Running

Start by supplying `{version}/{dataset}/rat_ids.txt` and `fastq_map.txt` for every selected **source dataset**. For existing datasets, copy these from `samples/rat_ids_current/rat_ids-{dataset}.txt` and `samples/fastq_map_current/fastq_map-{dataset}.txt`. The current list/map files are the source of the derived sample tracking tables, not outputs generated from those tables. New datasets require list/map inputs derived from their sequencing metadata. For a tissue containing a single source dataset, those files also serve as the tissue's files.

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

### Running one tissue or dataset

Keep the full `tissues` configuration and request specific output paths from the RatGTEx repository root. Commenting out tissues also changes the sample union used to build the shared genotype panel, which can trigger a panel rebuild and downstream processing. The stage runner selects all configured tissues; use Snakemake directly to select individual outputs.

For example, test alignment and identity QC on BLA:

```shell
snakemake v5/BLA/qc/BLA.star_stats.tsv --profile slurm -j10 -n
snakemake v5/BLA/qc/rna_to_geno_summary.tsv --profile slurm -j10 -n
```

Run the alignment command first, removing `-n` after reviewing the dry run. It produces coordinate-sorted BAMs, BAM indexes, and STAR statistics. Then run identity QC and review the results before proceeding to Pantry. Missing individual-sample VCFs are extracted from the shared panel as needed.

Alignment and identity QC use **source dataset names**: for NAcc, request NAcc1, NAcc2, NAcc3, and NAcc4 separately. Genotype hardcalls and downstream tissue analyses use **analysis tissue names**, such as NAcc.

Other individual output targets include:

```shell
# Tissue BED/BIM/FAM genotype files
snakemake v5/Brain/geno.bed --profile slurm -j10 -n

# Cross-modality phenotypes after Pantry phenotyping
snakemake v5/Brain/phenos/output/cross_modality.bed.gz --profile slurm -j10 -n

# One supplementary QTL output after Pheast QTL mapping
snakemake v5/Brain/Brain.expression.cis_qtl_signif.txt.gz --profile slurm -j10 -n
```

The last command requests one output from `qtl_extra`; the aggregate stage also generates expression p-values, aFC, and trans-QTL outputs. These targets consume existing Pantry products and do not launch Pantry. Run Pantry phenotyping and Pheast for individual tissues in their project directories, as described below.

Snakemake builds missing or outdated dependencies for each requested output. To resume processing all configured datasets later, use `python scripts/run_stage.py align` (or the corresponding stage). Without a target, Snakemake uses `default_stage`, initially `align`.

### Identity review before Pantry setup

Review identity matches manually. For unresolved samples, list their IDs plus a positive control in `{version}/{dataset}/qc/samples_without_matches.txt` and request `{version}/{dataset}/qc/all_rats_summary.tsv` with Snakemake. This builds the optional all-rats exon VCF as needed.

Record relabeling/removal decisions, update the current rat-ID lists and FASTQ maps and their release-local copies, then regenerate the derived sample tracking tables. Run `setup` to refresh combined-tissue files, then `genotypes` to refresh genotype inputs. Delete BAMs for changed identities and regenerate alignments and dependent outputs. Genotype-aware alignment makes rebuilding these BAMs necessary. Repeat identity QC before `pantry_setup`.

The unexpected-identity search uses an exon-only VCF containing all original round 11.3 samples. To prepare that shared reference separately:

```shell
snakemake geno/ratgtex_v5_round11_3.all_rats_exons.vcf.gz \
    geno/ratgtex_v5_round11_3.all_rats_exons.snps.txt --profile slurm -j10 -n
```

### Pantry projects and tissue assembly

`pantry_setup` creates `{version}/{name}/phenos` and `pheast` projects for selected source datasets and analysis tissues. It refuses existing Pantry project directories. It records the Pantry revision, working-tree status, and main configuration in `ratgtex_setup.json`.

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

Combined tissues reuse source-dataset sample-level products through symlinks, then compute tissue-dependent products in their own projects. Their maps contain absolute FASTQ paths, so mixed read lengths/layouts retain their original source quantification. Overlapping animal IDs across constituent datasets raise an error rather than silently deduplicating libraries. Decide how those libraries should be combined before proceeding.

### Workflow checks

```shell
python -m unittest discover -s tests -v
```

These checks exercise setup isolation, sample identity handling, config generation, and stage dispatch. Use Snakemake dry runs to validate the actual selected data and dependencies.

## Help

Email me or file an issue on this GitHub repository if anything isn't working, and I'll be happy to fix it.
