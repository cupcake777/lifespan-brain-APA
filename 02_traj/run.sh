#!/usr/bin/env bash
# Lifespan trajectory models, validation tests and manuscript figures (module 02).
# Needs APA_ROOT exported (source config/paths.local.sh with set -a first).
set -euo pipefail
: "${APA_ROOT:?source config/paths.local.sh first}"

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1

wait_for_all() {
  local status=0 pid
  for pid in "$@"; do
    wait "$pid" || status=1
  done
  return "$status"
}
cd "$(dirname "${BASH_SOURCE[0]}")"

# --- primary lifespan models ---

Rscript 01_cov.R
Rscript 02_support.R
Rscript 02_gam.R
Rscript 02_prediction.R
Rscript 03_dynamic.R
APA_CLUSTER_BOOT_ITER=100 Rscript 04_cluster.R

# --- primary downstream interpretation ---

Rscript 05_func.R

# --- validation and external tests ---

if (( TOTAL_CPUS > 1 )); then
  VALIDATION_CPUS=$(( (TOTAL_CPUS + 1) / 2 ))
  (
    export SLURM_CPUS_PER_TASK="$VALIDATION_CPUS"
    Rscript 08_brainspan.R
  ) &
  BRAINSPAN_PID=$!
  (
    export SLURM_CPUS_PER_TASK="$VALIDATION_CPUS"
    Rscript 10_batch9_holdout.R
  ) &
  HOLDOUT_PID=$!
  if ! wait_for_all "$BRAINSPAN_PID" "$HOLDOUT_PID"; then
    exit 1
  fi
else
  Rscript 08_brainspan.R
  Rscript 10_batch9_holdout.R
fi

if (( TOTAL_CPUS > 1 )); then
  Rscript 11_brainspan_heat.R &
  HEATMAP_PID=$!
  Rscript 12_ext_pattern.R &
  EXTERNAL_PID=$!
  Rscript fig_12_c2_scatter.R &
  DIAGNOSTIC_PID=$!
  if ! wait_for_all "$HEATMAP_PID" "$EXTERNAL_PID" "$DIAGNOSTIC_PID"; then
    exit 1
  fi
else
  Rscript 11_brainspan_heat.R
  Rscript 12_ext_pattern.R
  Rscript fig_12_c2_scatter.R
fi

# --- manuscript figures ---

FIGURE_JOBS="$TOTAL_CPUS"
figure_pids=()
for script in \
  fig_01b.R \
  fig_02.R \
  fig_03.R \
  fig_04.R \
  fig_05.R \
  fig_06.R \
  fig_07.R \
  fig_08.R \
  fig_09.R \
  fig_08_brainspan.R
do
  Rscript "$script" &
  figure_pids+=("$!")
  if (( ${#figure_pids[@]} >= FIGURE_JOBS )); then
    if ! wait_for_all "${figure_pids[@]}"; then
      exit 1
    fi
    figure_pids=()
  fi
done
if (( ${#figure_pids[@]} > 0 )); then
  wait_for_all "${figure_pids[@]}"
fi
