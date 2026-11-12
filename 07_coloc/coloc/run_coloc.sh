#!/usr/bin/env bash
# Run the EasyColoc managed pipeline; results go to data/07_coloc/coloc_run.
# Submit from $APA_ROOT: sbatch release/07_coloc/coloc/run_coloc.sh (create logs/ first).
#SBATCH --job-name=apa_coloc
#SBATCH --output=logs/coloc_%j.out
#SBATCH --error=logs/coloc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=48:00:00
#SBATCH --partition=CU

set -eu
: "${APA_ROOT:?set APA_ROOT (see release/config/paths.sh)}"

# --- config ---
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EASYCOLOC_ROOT="${EASYCOLOC_ROOT:-${HERE}/../easycoloc}"
RESULTS_DIR="${RESULTS_DIR:-${APA_ROOT}/data/07_coloc/coloc_run}"
GLOBAL_CFG="${GLOBAL_CFG:-${EASYCOLOC_ROOT}/config/global.yml}"
GWAS_CFG="${GWAS_CFG:-${EASYCOLOC_ROOT}/config/gwas.yml}"
QTL_CFG="${QTL_CFG:-${EASYCOLOC_ROOT}/config/qtl.yml}"
MAMBA_ENV_NAME="${MAMBA_ENV_NAME:-easycoloc}"
MICROMAMBA_BIN="${MICROMAMBA_BIN:-micromamba}"
N_CORES="${SLURM_CPUS_PER_TASK:-8}"

# --- env ---
export OMP_NUM_THREADS="${N_CORES}"
export MKL_NUM_THREADS="${N_CORES}"
export OPENBLAS_NUM_THREADS="${N_CORES}"
export EASYCOLOC_ASSIGN_RSID="always"

# --- run ---
mkdir -p "${RESULTS_DIR}" "${APA_ROOT}/logs"
cd "${EASYCOLOC_ROOT}"
"${MICROMAMBA_BIN}" run -n "${MAMBA_ENV_NAME}" \
  ./easycoloc run \
    --managed \
    --global "${GLOBAL_CFG}" \
    --gwas "${GWAS_CFG}" \
    --qtl "${QTL_CFG}" \
    --output-dir "${RESULTS_DIR}"
