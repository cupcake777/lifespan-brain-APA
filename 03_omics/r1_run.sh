#!/usr/bin/env bash
#SBATCH --job-name=figdata_r1
#SBATCH --partition=cu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/r1_run.%j.out
#SBATCH --error=logs/r1_run.%j.err

set -euo pipefail
Rscript "${APA_ROOT}/release/03_omics/r1_utr_expression.R"
