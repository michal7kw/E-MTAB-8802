#!/bin/bash
#SBATCH --job-name=analyze_rnaseq
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/analyze_rnaseq.err"
#SBATCH --output="logs/analyze_rnaseq.out"

# Stop on error
set -e
set -u
set -o pipefail

# Activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate snakemake

# Define working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/E-MTAB-8802"
cd $WORKDIR || { log_message "ERROR: Failed to change to working directory"; exit 1; }

Rscript analyze_rnaseq.R