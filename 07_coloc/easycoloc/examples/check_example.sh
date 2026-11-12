#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

bash examples/checks/check_cli.sh
bash examples/checks/check_report_web_cli.sh
Rscript examples/synthetic_plot_demo.R

echo "[EXAMPLE] all example checks passed"
