# RatGTEx Pipeline

This is the code used to process data for the [RatGTEx Portal](https://ratgtex.org). It is built on [Snakemake](https://snakemake.github.io/), a Python-based framework for reproducible data analysis. Since things might not work the first time with a new dataset, it makes it easy to process data iteratively. You can run part of the pipeline on a subset of the data, and then run it on the full dataset without regenerating existing files.

To add a new tissue to RatGTEx:
1. Make a directory within this one named as the tissue abbreviation.
2. Set up the software as described below.
3. Set up the raw data as described below.
4. Run snakemake as described below.

The main steps of the pipeline are:
1. Align RNA-Seq reads using [STAR](https://github.com/alexdobin/STAR).
2. Quantify gene expression using [RSEM](https://deweylab.github.io/RSEM/).
3. Map cis-eQTLs and trans-eQTLs using [tensorQTL](https://github.com/broadinstitute/tensorqtl) in various modes.
4. Calculate cis-eQTL effect size (allelic fold change) using [aFC.py](https://github.com/secastel/aFC).

## Setup

### Conda environment

Install these conda packages (I'm probably missing some):

- snakemake
- star
- rsem
- gtfparse
- tabix
- bioconductor-variantannotation
- plink2=2

For tensorQTL:
- tensorqtl (pip)
- rpy2
- bioconductor-qvalue
- fastparquet

For aFC:
- pysam
- statsmodels
- scikits-bootstrap

### Software

tensorQTL is available on pip, but install from GitHub for the latest version. I think this works:

```
cd tools
git clone git@github.com:broadinstitute/tensorqtl.git
pip install -e tensorqtl
```

To get `aFC.py`:

```
cd tools
git clone git@github.com:secastel/aFC.git
```

### Snakemake profile

When you run snakemake, you specify a profile that determines how steps get run. Here is the config file I use on a computing cluster with slurm scheduling:

`~/.config/snakemake/slurm/config.yaml`:

```yaml
use-conda: true
cluster: "sbatch -t {resources.walltime}:00:00 --mem={resources.mem_mb} -c {resources.cpus} {resources.partition} --mail-type=FAIL --mail-user=dmunro@scripps.edu"
default-resources: [walltime=1, mem_mb=4000, cpus=1, partition=""]
# partition should be e.g. "--partition=gpu"
```

You will then find resources specified for some of the snakemake rules, which are plugged into this command and automatically submitted as cluster jobs.

## Input files

### `geno/ratgtex.vcf.gz`

This is the merge of genotypes from all rats across tissues. This is done to unify the format of original genotypes, and because the portal expects a single set of SNPs, even if some tissues don't have genotypes for some SNPs.

### `geno/founders.vcf.gz`

These are the genotypes for the eight HS rat founder strains. They are used to calculate genotype covariates.

### `geno/sim_to_founders.txt`

These similarities to each founder strains are used as covariates. If this file exists already and new rats are added to `ratgtex.vcf.gz` for the new tissue, this file must be regenerated using `Snakefile`.

### FASTQ files

The FASTQ files can be in any accessible location. Currently only single-read RNA-Seq is implemented, but paired-end will also be supported soon.

### `{tissue}/fastq_map.txt`

A tab-delimited file with no header containing FASTQ file names and the rat IDs they correspond to. If multiple files map to the same ID, reads from all those files will be aligned into one BAM file. The file names can be absolute or relative paths in relation to this directory, and/or `Snakefile` can be edited to specify the path to the files.

### `{tissue}/rat_ids.txt`

A file listing the rat IDs for the dataset, one per line.

### Quality Control

Note that this pipeline currently does not include QC steps like filtering bad samples and correcting sample ID mixups. Instead, if QC has been performed previously e.g. in the dataset's original study, that information is reflected in `fastq_map.txt` as omitted or corrected rat IDs.

## Running

If everything is set up correctly, you can do a dry run:

`snakemake --profile slurm -j8 Eye/Eye.aFC.txt -n`

and if the steps seem to be what you expect, run it without the `-n` tag.

You may want to run the heavy raw data processing steps first, then move on once those are done. E.g. add the list of BAM files to the first rule in `Snakefile` and generate them:

`snakemake --profile slurm -j8`

## Help

There are likely a number of issues remaining with this pipeline, so file an issue on this GitHub repository if anything isn't working, and I'll be happy to fix it.
