#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

show_help() {
  cat <<'EOF'
Usage: tools/run_pipeline_managed.sh [--global PATH] [--gwas PATH] [--qtl PATH] [--output-dir PATH] [--shard-index N] [--shard-count N] [--shard-by locus|gwas]

Starts EasyColoc with managed logging, runtime heartbeat, manifest build,
monitor snapshot and completion check.
EOF
}

global_cfg=""
gwas_cfg=""
qtl_cfg=""
output_dir=""
shard_index=""
shard_count=""
shard_by=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --global)
      global_cfg="${2:?missing value for --global}"
      shift 2
      ;;
    --gwas)
      gwas_cfg="${2:?missing value for --gwas}"
      shift 2
      ;;
    --qtl)
      qtl_cfg="${2:?missing value for --qtl}"
      shift 2
      ;;
    --output-dir)
      output_dir="${2:?missing value for --output-dir}"
      shift 2
      ;;
    --shard-index)
      shard_index="${2:?missing value for --shard-index}"
      shift 2
      ;;
    --shard-count)
      shard_count="${2:?missing value for --shard-count}"
      shift 2
      ;;
    --shard-by)
      shard_by="${2:?missing value for --shard-by}"
      shift 2
      ;;
    --help|-h)
      show_help
      exit 0
      ;;
    *)
      echo "[RUN] unknown argument: $1" >&2
      show_help >&2
      exit 1
      ;;
  esac
done

if [[ -n "${global_cfg}" ]]; then
  export EASYCOLOC_GLOBAL_CONFIG="${global_cfg}"
fi
if [[ -n "${gwas_cfg}" ]]; then
  export EASYCOLOC_GWAS_CONFIG="${gwas_cfg}"
fi
if [[ -n "${qtl_cfg}" ]]; then
  export EASYCOLOC_QTL_CONFIG="${qtl_cfg}"
fi
if [[ -n "${shard_index}" ]]; then
  export EASYCOLOC_SHARD_INDEX="${shard_index}"
  export EASYCOLOC_SHARD_MODE="shard"
fi
if [[ -n "${shard_count}" ]]; then
  export EASYCOLOC_SHARD_COUNT="${shard_count}"
  export EASYCOLOC_SHARD_MODE="shard"
fi
if [[ -n "${shard_by}" ]]; then
  export EASYCOLOC_SHARD_BY="${shard_by}"
fi

if [[ -z "${output_dir}" ]]; then
  output_dir="$(Rscript -e 'source("src/utils_config.R"); cfg <- easycoloc_read_configs(); cat(normalizePath(cfg$global$output_dir, mustWork = FALSE))')"
else
  export EASYCOLOC_OUTPUT_DIR="${output_dir}"
fi

mkdir -p "${output_dir}"

timestamp="$(date +%Y%m%d_%H%M%S)"
if [[ "${EASYCOLOC_SHARD_MODE:-single}" == "shard" ]]; then
  shard_index_label="${EASYCOLOC_SHARD_INDEX:-1}"
  shard_count_label="${EASYCOLOC_SHARD_COUNT:-1}"
  shard_by_label="${EASYCOLOC_SHARD_BY:-locus}"
  log_file="${output_dir}/run_easycoloc_shard_${shard_by_label}_${shard_index_label}of${shard_count_label}_${timestamp}.log"
  run_label="managed_shard_${shard_by_label}_${shard_index_label}of${shard_count_label}_${timestamp}"
else
  log_file="${output_dir}/run_easycoloc_full_${timestamp}.log"
  run_label="managed_${timestamp}"
fi

total_gwas="$(Rscript -e 'source("src/utils_config.R"); cfg <- easycoloc_read_configs(); cat(length(cfg$gwas$datasets))')"

echo "[RUN] log file: ${log_file}"
echo "[RUN] starting EasyColoc pipeline"
if [[ "${EASYCOLOC_SHARD_MODE:-single}" == "shard" ]]; then
  echo "[RUN] shard: ${EASYCOLOC_SHARD_BY:-locus} ${EASYCOLOC_SHARD_INDEX:-1}/${EASYCOLOC_SHARD_COUNT:-1}"
fi

EASYCOLOC_LOG_FILE="${log_file}" \
EASYCOLOC_RUN_LABEL="${run_label}" \
EASYCOLOC_PARENT_PID="$$" \
Rscript run_coloc.R > "${log_file}" 2>&1 &
run_pid=$!

cleanup() {
  if kill -0 "${run_pid}" 2>/dev/null; then
    kill "${run_pid}" 2>/dev/null || true
  fi
}
trap cleanup INT TERM

render_progress() {
  local processed_gwas current_gwas abf_count rds_count plot_count filled bar_width
  processed_gwas="$(grep -c '^Processing GWAS Dataset:' "${log_file}" 2>/dev/null || true)"
  current_gwas="$(grep '^Processing GWAS Dataset:' "${log_file}" 2>/dev/null | tail -n 1 | sed 's/^Processing GWAS Dataset: //' || true)"
  abf_count="$(if [[ -d "${output_dir}/abf" ]]; then find "${output_dir}/abf" -maxdepth 1 -type f -name '*_locus_results.csv' | wc -l | tr -d ' '; else printf '0'; fi)"
  rds_count="$(if [[ -d "${output_dir}/rds" ]]; then find "${output_dir}/rds" -maxdepth 1 -type f -name '*.rds' | wc -l | tr -d ' '; else printf '0'; fi)"
  plot_count="$(if [[ -d "${output_dir}/plots" ]]; then find "${output_dir}/plots" -maxdepth 1 -type f \( -name '*.pdf' -o -name '*.png' \) | wc -l | tr -d ' '; else printf '0'; fi)"
  bar_width=24
  if [[ "${total_gwas}" -gt 0 ]]; then
    filled=$(( processed_gwas * bar_width / total_gwas ))
  else
    filled=0
  fi
  printf -v bar '%*s' "${filled}" ''
  bar="${bar// /#}"
  printf -v pad '%*s' $(( bar_width - filled )) ''
  pad="${pad// /-}"
  printf '\r[%-24s] GWAS %s/%s | current: %s | abf:%s rds:%s plots:%s' \
    "${bar}${pad}" "${processed_gwas}" "${total_gwas}" "${current_gwas:-starting}" "${abf_count}" "${rds_count}" "${plot_count}"
}

while kill -0 "${run_pid}" 2>/dev/null; do
  render_progress
  sleep 5
done

set +e
wait "${run_pid}"
exit_code=$?
set -e
printf '\n'

echo "[RUN] pipeline exited with code ${exit_code}"
if [[ "${EASYCOLOC_SHARD_MODE:-single}" == "shard" ]]; then
  echo "[RUN] shard mode detected; skipping global manifest and completion check"
  exit "${exit_code}"
fi
echo "[RUN] building output manifest"
if Rscript tools/build_output_manifest.R "${output_dir}"; then
  echo "[RUN] manifest written"
else
  echo "[RUN] manifest build failed"
fi

echo "[RUN] collecting runtime monitor snapshot"
if Rscript tools/monitor_easycoloc.R "${output_dir}"; then
  echo "[RUN] monitor snapshot written"
else
  echo "[RUN] monitor snapshot failed"
fi

echo "[RUN] checking completion"
if Rscript tools/check_run_completion.R "${output_dir}"; then
  echo "[RUN] completion check passed"
  completion_exit_code=0
else
  echo "[RUN] completion check failed"
  completion_exit_code=1
fi

if [[ "${completion_exit_code}" -eq 0 ]]; then
  if [[ "${exit_code}" -ne 0 ]]; then
    echo "[RUN] pipeline process exited non-zero, but completion check passed; treating run as complete"
  fi
  exit 0
fi

if [[ "${exit_code}" -ne 0 ]]; then
  exit "${exit_code}"
fi
exit "${completion_exit_code}"
