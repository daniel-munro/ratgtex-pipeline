# RatGTEx Pipeline

This is the code used to process data for the [RatGTEx Portal](https://ratgtex.org). It is built on [Snakemake](https://snakemake.github.io/), a Python-based framework for reproducible data analysis. Since things might not work the first time with a new dataset, it makes it easy to process data iteratively. You can run part of the pipeline on a subset of the data, and then run it on the full dataset without regenerating existing files.

To add a new tissue to RatGTEx:
1. Make a directory within this one named as the tissue abbreviation.
2. Set up the software as described below.
3. Set up the raw data as described below.
4. Run snakemake as described below.

The main steps of the pipeline are:
1. Align RNA-Seq reads using [STAR](https://github.com/alexdobin/STAR).
2. Checking for and fixing sample mixups.
3. Quantify gene expression using [RSEM](https://deweylab.github.io/RSEM/).
4. Map cis-eQTLs and trans-eQTLs using [tensorQTL](https://github.com/broadinstitute/tensorqtl) in various modes.
5. Calculate cis-eQTL effect size (allelic fold change) using [aFC.py](https://github.com/secastel/aFC).
6. Map cis-sQTLs and trans-sQTLs using [regtools](https://regtools.readthedocs.io/en/latest/), [leafCutter](http://davidaknowles.github.io/leafcutter/), and tensorQTL.

Snakemake automatically links the pipeline together based on input and output files. Here is how all the steps link together:

![RatGTEx pipeline rulegraph](test/rulegraph.png)

Here's a more detailed version showing the inputs and outputs:

![RatGTEx pipeline filegraph](test/filegraph.png)

## Setup

### Conda environment

Install conda and [add the bioconda channel](https://bioconda.github.io/user/install.html#set-up-channels).

Install these packages using conda, e.g. `conda install snakemake`, or pip if indicated in parentheses:

- snakemake
- star
- rsem
- gtfparse
- tabix
- bioconductor-variantannotation
- bioconductor-impute
- plink2=2
- pandas
- bx-python
- scipy

For QC:
- bcftools
- gatk4

For tensorQTL:
- tensorqtl (pip)
- rpy2
- bioconductor-qvalue
- fastparquet

For aFC:
- pysam
- statsmodels
- scikits-bootstrap

For splice QTLs:
- regtools
- scikit-learn
- r-argparser

(There's probably more that I've forgotten.)

### Other software

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
default-resources: [walltime=4, mem_mb=4000, cpus=1, partition=""]
# partition should either be empty string for default or something like "--partition=gpu"
latency-wait: 60
cluster-cancel: scancel
```

Resources are specified for some of the snakemake rules, which are plugged into this command and automatically submitted as cluster jobs.

On TSCC, which uses TORQUE scheduling, the jobs run in the default conda environment rather than the one that is active when executing snakemake. So if you want to use an environment other than `base`, e.g. `ratgtex`, add `conda activate ratgtex` to the end of your `~/.bashrc` file. This is just a workaround since it affects which environment loads on login too, and must be modified if you want to run jobs for another project. Let me know if you find a way to either specify the environment in the Snakefile or snakemake command, or always run jobs using the active environment.

### Pan-tissue input files

I'll provide reference files too big for this repo in `tscc:/home/dmunro/ratgtex`, which you can copy or make symbolic links to.

#### `ref_rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel.{dict,fa,fa.fai}` and `ref_rn7/Rattus_norvegicus.mRatBN7.2.dna.toplevel.{dict,fa,fa.fai}`

Rat genome files from Ensembl.

#### `ref_rn6/Rattus_norvegicus.Rnor_6.0.99.gtf` and `ref_rn7/Rattus_norvegicus.mRatBN7.2.108.gtf`

Gene annotations from Ensembl.

#### `geno_rn{6/7}/all_rats_exons.vcf.gz`

This contains the genotypes for exon regions from the thousands of HS rats the Palmer Lab has collected. It is used to try to find a match for any RNA-seq samples that fail the sample mixup QC step and don't match any of the tissue's original cohort genotypes. It can be generated with `scripts/genotypes_rn{6/7}.sh`.

#### `config.yaml`

This configuration file contains parameters about the datasets and is used by snakemake. Parameters for each tissue are grouped under tissue names, e.g.:

```yaml
# IL, LHb, NAcc, OFC, and PL datsets
IL:
  read_length: 100
  fastq_path: "../gpfs/fastq"
  paired_end: false
  geno_dataset: "IL_LHb_NAcc_OFC_PL"
LHb:
  read_length: 100
  fastq_path: "../gpfs/fastq"
  paired_end: false
  geno_dataset: "IL_LHb_NAcc_OFC_PL"
...
```

Some of these are inherent to the data, while others, e.g. `fastq_path`, may need to be edited to point to the correct location.

### Dataset-specific input files

#### FASTQ files

The FASTQ files can be in any accessible location. Currently only single-read RNA-Seq is implemented, but paired-end will also be supported soon.

#### `rn{6/7}/{tissue}/fastq_map.txt`

A tab-delimited file with no header containing the paths to each FASTQ file and the rat IDs they correspond to. Or, for paired-end reads, each row contains the first FASTQ path, second FASTQ path, and rat ID per file pair.
- If multiple files map to the same ID, i.e. the ID appears in multiple rows, reads from those files will be aligned into one BAM file.
- You can use the `fastq_path` parameter in `config.yaml` to specify the encompassing directory as an absolute or relative path. That way `fastq_map.txt` can just contain the remainder of the path to each file (including any subdirectories as necessary).
- Any listed files whose rat IDs are not in `rn{6/7}/{tissue}/rat_ids.txt` will be ignored.

#### `rn{6/7}/{tissue}/rat_ids.txt`

A file listing the rat IDs for the dataset, one per line. This list determines which samples are included in the processing.

#### `geno_rn{6/7}/{dataset}.vcf.gz`

A VCF file containing the genotypes for one or more tissues. If multiple tissues came from the same project and have overlapping sets of individuals, they use the same VCF file. These are created using `scripts/genotypes_rn{6/7}.sh`, which ensures that REF alleles match the reference genome and that the VCFs are otherwise compatible with the pipeline. Specify the dataset name in a tissue- or dataset-specific config file as described below.

## Running

Edit `config.yaml` in this directory so that the tissue(s) you want to process are present and have correct parameters. Unlike the Snakemake config file, which specifies how jobs are run, this one contains parameters for the tissues/datasets such as read length and directory where FASTQ files can be found.

### QC

#### Pre-run checks

Before running Snakemake, run `python3 scripts/qc_init_check.py {rn6/rn7} {tissuename}`, which checks the input data and config for issues.

#### Sample mixup checks

The way to do sample mixup testing is to generate the mixup checking outputs using Snakemake, which will generate the BAM files as dependencies if needed. Examine the outputs to identify samples that need to be relabeled (e.g. if two labels get swapped) or removed.

- To relabel a sample, edit the ID in the 2nd column of `fastq_map.txt` for all of its FASTQ files so that its BAM file gets labeled correctly. You'll then need to regenerate the BAM file since it will now use the correct VCF individual as input to STAR.
- To remove a sample, remove its ID from `rat_ids.txt` and delete its BAM and any other generated files.

Before removing samples, run the second stage of sample mixup checking, which tests the RNA-seq samples that still don't have matches against 6000+ rat genotypes to see if a match can be found. To do this, list the mismatched samples in `rn{6/7}/{tissue}/qc/samples_without_matches.txt`, along with an OK sample as a positive control (if that sample is included in the all-rat VCF). Then generate `rn{6/7}/{tissue}/qc/all_rats_summary.tsv` and use any additional matches found. This will probably require adding the new matching genotypes to the VCF file (see `scripts/genotypes_rn{6/7}.sh`).

### Continue

After these corrections, continue with the pipeline.

You may want to run a subset of the heavy raw data processing steps first, then move on once those are done. E.g. add the first 10 BAM files to the first rule (called 'all') in `Snakefile` and generate them:

`snakemake --profile slurm -j10`

Use the `-n` dry run tag to make sure things seem to be set up correctly before running.

## Help

There are likely a number of issues remaining with this pipeline, so email me or file an issue on this GitHub repository if anything isn't working, and I'll be happy to fix it.
