#!/usr/bin/env bash
# Run the EasyColoc managed pipeline; output goes to one run folder under data/07_coloc (set RUN_NAME).
# Submit after sourcing config/paths.local.sh.
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
: "${APA_ROOT:?set APA_ROOT (see config/paths.local.sh)}"
: "${EASYCOLOC_ROOT:?set EASYCOLOC_ROOT to an external EasyColoc installation}"

# --- config ---
RUN_DIR="${APA_ROOT}/data/07_coloc/${RUN_NAME:?set RUN_NAME (new run folder name)}"
GLOBAL_CFG="${EASYCOLOC_ROOT}/config/global.yml"
GWAS_CFG="${EASYCOLOC_ROOT}/config/gwas.yml"
QTL_CFG="${EASYCOLOC_ROOT}/config/qtl.yml"
MAMBA_ENV_NAME="${MAMBA_ENV_NAME:-easycoloc}"
MICROMAMBA_BIN="${MICROMAMBA_BIN:-micromamba}"
N_CORES="${SLURM_CPUS_PER_TASK:-8}"

# --- env ---
export OMP_NUM_THREADS="${N_CORES}"
export MKL_NUM_THREADS="${N_CORES}"
export OPENBLAS_NUM_THREADS="${N_CORES}"
export EASYCOLOC_ASSIGN_RSID="always"

# --- run ---
mkdir -p "${RUN_DIR}" "${APA_ROOT}/logs"
cd "${EASYCOLOC_ROOT}"
"${MICROMAMBA_BIN}" run -n "${MAMBA_ENV_NAME}" \
  ./easycoloc run \
    --managed \
    --global "${GLOBAL_CFG}" \
    --gwas "${GWAS_CFG}" \
    --qtl "${QTL_CFG}" \
    --output-dir "${RUN_DIR}"
