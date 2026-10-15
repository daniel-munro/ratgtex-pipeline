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

RatGTEx orchestration requires Python with PyYAML and Snakemake with the Slurm executor plugin. Alignment/QC and supplementary QTL rules also use the tools invoked in `steps/` (including STAR, samtools, bcftools, PLINK, GATK, tensorQTL, pandas, NumPy, fastparquet, bgzip, and tabix). Use your installed scientific software environment for those rules. Pantry supplies its own environment specifications; see its documentation for phenotyping and Pheast dependencies. There is no root RatGTEx `environment.yml`.

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

#### `config.yaml`

`config.yaml` defines source **datasets**, analysis **tissues**, genotype panels, modalities, reference files, and Pantry settings.

- A **dataset** is a tissue-specific sequencing collection from a source study, such as NAcc1.
- A **tissue** is the tissue-level analysis assembled from one or more datasets, such as NAcc from NAcc1–4 or Brain from Brain alone.
- A **genotype panel** identifies the VCF used for a dataset, through `genotype_panel: <name>` and `geno/<name>.vcf.gz`. This is independent of RNA-seq dataset names.
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
    genotype_panel: "your_genotype_panel"
  # Define each NAcc source dataset here too.
```

Only datasets referenced by `tissues` are selected. Tissue membership replaces the old `run`/`merged_tissues` lists and numeric-suffix convention. `modalities` defines each phenotype's Pantry producer and grouping; `cross_modality: true` combines all configured modalities for grouped cis-QTL analysis. TWAS uses the individual modalities because current Pheast restricts cross-modality analysis to QTL mapping.

`pantry.source` identifies the Pantry checkout copied during setup. `pantry.phenotyping` and `pantry.pheast` supply settings for generated project configs; `pantry.merge_products` declares the reusable sample-level products for each producer. Optional producer settings can be supplied under `pantry.phenotyping.modality_groups.<producer>`; output file lists are generated from `modalities`. Parameters still inside Pantry rules are not moved by this refactor.

The checked-in configuration still describes v4 inputs. Choose the v5 version, genotype panel, datasets/tissues, and modalities before setting up a new release. Paths in the main config are relative to the RatGTEx repository; generated Pantry configs use absolute input paths. Restore or configure the reference FASTA/GTF and TWAS SNP list before running new processing; the paths currently configured for those inputs were absent during this review.

### Dataset-specific input files

#### FASTQ files

The RNA-seq FASTQ files can be single-read or paired-end. There can be more than one file or file pair per sample, and all will be used.

#### `{version}/{dataset}/fastq_map.txt`

A tab-delimited file with no header containing the paths to each FASTQ file and the rat IDs they correspond to. Or, for paired-end reads, each row contains the first FASTQ path, second FASTQ path, and rat ID per file pair.

- If multiple files map to the same ID, i.e. the ID appears in multiple rows, reads from those files will be aligned into one BAM file.
- You can use the `fastq_path` parameter in `config.yaml` to specify the encompassing directory as an absolute or relative path. That way `fastq_map.txt` can just contain the remainder of the path to each file (including any subdirectories as necessary).
- Any listed files whose rat IDs are not in `{version}/{dataset}/rat_ids.txt` will be ignored.
- If a dataset includes a mix of single-end and paired-end samples, omit `paired_end` from the dataset config so the pipeline infers pairing per row.

#### `{version}/{dataset}/rat_ids.txt`

A file listing the rat IDs for the dataset, one per line. This list determines which samples are included in the processing.

#### `geno/{genotype_panel}.vcf.gz`

A VCF file containing the genotypes for one or more tissues. If multiple tissues came from the same project and have overlapping sets of individuals, they use the same VCF file. These are created using `scripts/setup/genotypes_{version}.sh`, which ensures that REF alleles match the reference genome and that the VCFs are otherwise compatible with the pipeline. Specify the panel name in each source dataset's `genotype_panel` setting.

## Running

From the RatGTEx repository, use the stage runner for setup and the RatGTEx processing steps:

```shell
python scripts/run_stage.py align --config config.yaml -n
```

`--profile slurm` and `--jobs 10` are defaults. Remove `-n` to execute a processing stage. `--plan` prints commands without invoking Snakemake or changing files; `setup -n` and `configure -n` preview their actions without writing files. A processing dry run requires its input data and project setup to be present. Pantry phenotyping and Pheast are run manually in individual project directories; the stage runner never launches them.

| Stage            | Action                                                                                                                                        |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `align`          | Align selected source datasets and produce BAM indexes and STAR summaries.                                                                    |
| `identity_qc`    | Produce RNA/genotype identity summaries, generating alignments as needed.                                                                     |
| `setup`          | Copy new self-contained Pantry projects, generate configs, assemble tissue sample lists/maps, and link reusable inputs.                       |
| `configure`      | Regenerate phenotyping and named QTL/TWAS configs after main-config parameter changes; preserve previous configs in `config_history/`.        |
| `genotypes`      | Generate tissue PLINK files.                                                                                                                  |
| `cross_modality` | Rebuild combined phenotype tables from all configured modalities.                                                                             |
| `qtl_extra`      | Generate RatGTEx significant-pair, expression p-value, aFC, and trans-QTL outputs after Pheast QTL mapping. Requires the expression modality. |

For processing stages, the runner invokes the main RatGTEx Snakefile once with the selected profile. Setup and configuration can prepare all selected projects together, but each Pantry project has its own independent Snakemake invocation, logs, and scheduler. Launch several projects in separate sessions when desired, and troubleshoot or rerun an individual project directly.

The RatGTEx-native targets are also available directly:

```shell
snakemake align --profile slurm -j10 -n
snakemake identity_qc --profile slurm -j10 -n
snakemake genotypes --profile slurm -j10 -n
snakemake cross_modality --profile slurm -j10 -n
snakemake qtl_extra --profile slurm -j10 -n
```

Without a target, Snakemake uses `default_stage` (initially `align`). Its supported values are `align`, `identity_qc`, `genotypes`, and `qtl_extra`. No changes to `rule all` are needed. To process a subset, use a config containing only the desired tissues, or request individual output paths directly with Snakemake.

### Identity review before phenotyping

Prepare source-dataset `rat_ids.txt` and `fastq_map.txt` inputs, then run `align` and `identity_qc`. Automatic generation of source-dataset sample files from tracking information remains a separate Todo.

Review identity matches manually. For unresolved samples, list their IDs plus a positive control in `{version}/{dataset}/qc/samples_without_matches.txt` and request `{version}/{dataset}/qc/all_rats_summary.tsv` with Snakemake. Update genotypes if newly identified animals must be added.

Record relabeling/removal decisions in the sample tracking information and update the processing maps/lists accordingly. Delete BAMs for changed identities and regenerate alignments and dependent outputs. Repeat QC before continuing. Genotype-aware alignment makes rebuilding these BAMs necessary.

After identity decisions are settled, run `setup`. The subsequent handoffs are explicit:

1. Run Pantry phenotyping manually for each source dataset. Independent datasets can run concurrently.
2. For a merged tissue, finish its constituent datasets' phenotyping first, then run the merged tissue's Pantry project manually. Unrelated tissues can proceed independently.
3. Run the RatGTEx `cross_modality` stage after tissue phenotyping, and `genotypes` before Pheast. Genotype preparation can run independently of phenotyping once tissue sample lists are ready.
4. Run Pheast QTL and/or TWAS manually for each tissue. QTL uses cross-modality outputs when enabled; TWAS uses the individual modalities.
5. Run RatGTEx `qtl_extra` after the relevant Pheast QTL results and covariates are ready.

The aggregate RatGTEx targets cover all configured tissues. To proceed with one ready tissue while others are still running, request its output files directly; for example:

```shell
snakemake v5/Brain/phenos/output/cross_modality.bed.gz --profile slurm -j10 -n
snakemake v5/Brain/geno.bed --profile slurm -j10 -n
snakemake v5/Brain/Brain.expression.cis_qtl_signif.txt.gz --profile slurm -j10 -n
```

Remove `-n` when ready to execute each command. These RatGTEx targets consume existing Pantry products and do not launch Pantry to generate missing products.

### Pantry projects and tissue assembly

`setup` creates `{version}/{name}/phenos` and `pheast` projects for selected source datasets and analysis tissues. It refuses existing Pantry project directories, protecting previous code and results; do not use it to migrate the existing v4 copies. It records the Pantry revision, working-tree status, and main configuration in `ratgtex_setup.json`.

Run phenotyping directly inside the project you want to process:

```shell
cd v5/Brain/phenos
snakemake --profile slurm -j20 -n
```

Remove `-n` to run. Start another tissue in a separate session to run it concurrently. For merged NAcc, run the same command inside `v5/NAcc/phenos` after NAcc1–4 have finished.

Each Pheast project has `config.qtl.yml`/`Snakefile.qtl` and `config.twas.yml`/`Snakefile.twas`. Select the desired entry point yourself. Initial `config.yml` defaults to QTL; subsequent `configure` calls leave that file alone so it can be used for additional interactive analyses. For example:

```shell
cd v5/Brain/pheast
snakemake --snakefile Snakefile.qtl --profile slurm -j10 -n
# Or, for TWAS:
snakemake --snakefile Snakefile.twas --profile slurm -j10 -n
```

Use `configure` for parameter changes after setup. Changes to tissue membership, sample lists, or merge-product links require reviewing/rebuilding the corresponding setup inputs and links; `configure` only regenerates configuration files, not tissue membership or symlinks.

Combined tissues reuse source-dataset sample-level products through symlinks, then compute tissue-dependent products in their own projects. Their maps contain absolute FASTQ paths, so mixed read lengths/layouts retain their original source quantification. Overlapping animal IDs across constituent datasets raise an error rather than silently deduplicating libraries. Decide how those libraries should be combined before proceeding. The old version-specific setup scripts remain historical documentation; the stage runner replaces them for new project setup.

### Workflow checks

```shell
python -m unittest discover -s tests -v
```

These checks exercise setup isolation, sample identity handling, config generation, and stage dispatch. Use Snakemake dry runs to validate the actual selected data and dependencies.

## Help

Email me or file an issue on this GitHub repository if anything isn't working, and I'll be happy to fix it.
