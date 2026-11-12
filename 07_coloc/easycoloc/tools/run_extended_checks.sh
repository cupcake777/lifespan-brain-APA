#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${repo_root}"

echo "[EXTENDED] Parsing core scripts"
Rscript tools/checks/check_parse.R

echo "[EXTENDED] Running PLINK clump smoke test"
Rscript tools/checks/smoke_test_plink_clump.R

echo "[EXTENDED] Running standard smoke checks"
./tools/run_smoke_checks.sh

echo "[EXTENDED] All extended checks passed"
