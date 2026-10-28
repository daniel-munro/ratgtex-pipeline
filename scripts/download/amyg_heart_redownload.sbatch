#!/bin/bash
#SBATCH --job-name=ena_Amyg_Heart_retry
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH --account=gdml
#SBATCH --partition=cpu-gdml
#SBATCH --output=fastq/Amyg_Heart/redownload-%j.log

set -euo pipefail

# Submit from the RatGTEx repository root. Run only one downloader at a time.
python -u scripts/download/amyg_heart.py redownload \
  --failure-log fastq/Amyg_Heart/redownload-8355123.log \
  --workers 4
