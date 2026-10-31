#!/usr/bin/env bash
set -euo pipefail

target="${1:-}"
case "$target" in
  figure_[1-6]|extended_data_[1-9]|main|extended_data|all) ;;
  *)
    printf 'Usage: %s {figure_1..figure_6|extended_data_1..extended_data_9|main|extended_data|all}\n' "$0" >&2
    exit 2
    ;;
esac

: "${APA_CODE_ROOT:?Set APA_CODE_ROOT (source config/paths.local.sh)}"
: "${APA_ROOT:?Set APA_ROOT (source config/paths.local.sh)}"
export APA_FIGURE_DIR="${APA_FIGURE_DIR:-${APA_CODE_ROOT}/outputs/figures}"
mkdir -p "${APA_FIGURE_DIR}"

run_r() {
  local script="$1"
  printf 'Rendering %s\n' "$script"
  Rscript "${APA_CODE_ROOT}/${script}"
}

run_python() {
  local script="$1"
  printf 'Rendering %s\n' "$script"
  python "${APA_CODE_ROOT}/${script}"
}

render_main() {
  local n="$1"
  case "$n" in
    1) run_r figure_generation/main_figures/figure_1.R ;;
    2) run_python figure_generation/panels/figure_2d_gene_model.py; run_r figure_generation/main_figures/figure_2.R ;;
    3) run_python figure_generation/panels/figure_3d_cep41_locus.py; run_r figure_generation/main_figures/figure_3.R ;;
    4) run_python figure_generation/panels/figure_4d_mechanism.py; run_r figure_generation/main_figures/figure_4.R ;;
    5) run_r figure_generation/main_figures/figure_5.R ;;
    6) run_r figure_generation/main_figures/figure_6.R ;;
    *) printf 'Unknown main figure: %s\n' "$n" >&2; exit 2 ;;
  esac
}

render_ed() {
  local n="$1"
  [[ "$n" =~ ^[1-9]$ ]] || { printf 'Unknown Extended Data figure: %s\n' "$n" >&2; exit 2; }
  run_r "figure_generation/extended_data_figures/extended_data_${n}.R"
}

case "$target" in
  figure_[1-6]) render_main "${target#figure_}" ;;
  extended_data_[1-9]) render_ed "${target#extended_data_}" ;;
  main) for n in {1..6}; do render_main "$n"; done ;;
  extended_data) for n in {1..9}; do render_ed "$n"; done ;;
  all)
    for n in {1..6}; do render_main "$n"; done
    for n in {1..9}; do render_ed "$n"; done
    ;;
esac
