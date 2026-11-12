#!/usr/bin/env bash
# Slurm wrapper for the lifespan trajectory run (module 02).
# Submit from APA_ROOT: sbatch release/02_traj/run.slurm
# Run `mkdir -p $APA_ROOT/logs` before sbatch; --output/--error do not create it.

#SBATCH --job-name=apa02
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/02_traj_%j.out
#SBATCH --error=logs/02_traj_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
set -a
source release/config/paths.local.sh
set +a
bash release/02_traj/run.sh
