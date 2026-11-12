#!/usr/bin/env bash
# Slurm wrapper for the protein validation run (module 03).
# Submit from APA_ROOT: sbatch release/03_omics/run.slurm

#SBATCH --job-name=apa03_protein
#SBATCH --partition=cu
#SBATCH --cpus-per-task=16
#SBATCH --mem=192G
#SBATCH --time=72:00:00
#SBATCH --output=logs/03_omics_%j.out
#SBATCH --error=logs/03_omics_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source release/config/paths.local.sh

Rscript release/03_omics/prep_matched.R
bash release/03_omics/run_protein.sh
