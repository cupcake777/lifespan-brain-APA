#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

run_step() {
  local label="$1"
  shift
  echo "[CHECK] ${label}"
  "$@"
}

run_step "Parsing core scripts" Rscript tools/checks/check_parse.R
run_step "Running plotting smoke check" Rscript tools/checks/smoke_test_plotting.R
run_step "Running summary/report smoke check" Rscript tools/checks/smoke_test_summary_report.R
run_step "Running runtime monitor smoke check" Rscript tools/checks/smoke_test_runtime_monitor.R
run_step "Running output-dir tool smoke check" Rscript tools/checks/smoke_test_output_dir_tools.R
run_step "Running doctor smoke check" Rscript tools/checks/smoke_test_doctor.R
run_step "Running reference listing smoke check" Rscript tools/checks/smoke_test_refs.R
run_step "Running reference bootstrap smoke check" Rscript tools/checks/smoke_test_bootstrap_refs.R
run_step "Running toy demo bootstrap smoke check" Rscript tools/checks/smoke_test_bootstrap_demo.R
run_step "Running 1KG setup smoke check" Rscript tools/checks/smoke_test_1kg_setup.R
run_step "Running GTEx bootstrap smoke check" Rscript tools/checks/smoke_test_gtex_bootstrap.R
run_step "Running LD cache smoke check" Rscript tools/checks/smoke_test_ld_cache.R
run_step "Running LD timeout smoke check" Rscript tools/checks/smoke_test_ld_timeout.R
run_step "Running chunked harmonization smoke check" Rscript tools/checks/smoke_test_chunked_harmonization.R
run_step "Running report web data smoke check" Rscript tools/checks/smoke_test_report_web_data.R
run_step "Running CLI example check" bash examples/checks/check_cli.sh
run_step "Running report web CLI example check" bash examples/checks/check_report_web_cli.sh
run_step "Running minimal synthetic demo" Rscript examples/synthetic_plot_demo.R

echo "[CHECK] All smoke checks passed"
