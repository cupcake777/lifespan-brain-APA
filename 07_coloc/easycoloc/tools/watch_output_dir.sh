#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

show_help() {
  cat <<'EOF'
Usage: tools/watch_output_dir.sh OUTPUT_DIR [INTERVAL_SECONDS] [LOG_FILE]

Poll an EasyColoc results directory repeatedly and append monitor snapshots to
a local log file.
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  show_help
  exit 0
fi

output_dir="${1:-}"
interval="${2:-60}"
log_file="${3:-}"

if [[ -z "${output_dir}" ]]; then
  echo "[WATCH] OUTPUT_DIR is required" >&2
  show_help >&2
  exit 1
fi

if ! [[ "${interval}" =~ ^[0-9]+$ ]] || [[ "${interval}" -lt 1 ]]; then
  echo "[WATCH] INTERVAL_SECONDS must be a positive integer" >&2
  exit 1
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
if [[ -z "${log_file}" ]]; then
  safe_name="$(basename "${output_dir}" | tr -c '[:alnum:]_.-' '_')"
  log_file="${repo_root}/logs/monitor/${safe_name}_watch_${timestamp}.log"
fi

mkdir -p "$(dirname "${log_file}")"

echo "[WATCH] output_dir: ${output_dir}"
echo "[WATCH] interval_seconds: ${interval}"
echo "[WATCH] log_file: ${log_file}"

while true; do
  {
    printf '===== %s =====\n' "$(date '+%Y-%m-%d %H:%M:%S')"
    Rscript tools/monitor_easycoloc.R --no-write-snapshot "${output_dir}"
    printf '\n'
  } >> "${log_file}" 2>&1
  sleep "${interval}"
done
