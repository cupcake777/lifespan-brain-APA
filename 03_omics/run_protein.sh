#!/usr/bin/env bash
# Protein validation steps for module 03: producers 01-10 and the figure steps.
# Writes to data/03_omics; reads the module 02 trajectories from data/02_traj.
: "${APA_ROOT:?source config/paths.local.sh}"

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_dir="${APA_ROOT}/data/03_omics"
# --- config ---

export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1

# --- helpers ---
wait_for_all() {
  local status=0 pid
  for pid in "$@"; do
    wait "${pid}" || status=1
  done
  return "${status}"
}

# --- producers ---
cd "${project_dir}"

Rscript "${script_dir}/01_universe.R"
Rscript "${script_dir}/02_fit_traj.R"
Rscript "${script_dir}/03_pairs.R"
Rscript "${script_dir}/03_integrate.R"

Rscript "${script_dir}/04_enrichment.R" &
pid_enrichment=$!
Rscript "${script_dir}/05_gam_sensitivity.R" &
pid_gam=$!
wait_for_all "${pid_enrichment}" "${pid_gam}"

Rscript "${script_dir}/10_gam_shape.R"

# --- figures ---
Rscript "${script_dir}/06_plot_cross.R" &
pid_cross_layer=$!
Rscript "${script_dir}/07_story_data.R" &
pid_story_data=$!
Rscript "${script_dir}/09_triad_overview.R" &
pid_triad=$!
wait_for_all "${pid_cross_layer}" "${pid_story_data}" "${pid_triad}"

Rscript "${script_dir}/08_plot_story.R"
