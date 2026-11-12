#!/usr/bin/env bash
# End-to-end synthetic tests for the module 02 chain: generator determinism, the runnable chain
# (run_e2e.sh), and the planted signal check. Exit 0 = all pass. Writes only to a temp dir.
# Env: APA_E2E_PYTHON (python with numpy, pandas), APA_E2E_RSCRIPT (Rscript with bluster for 04).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${APA_E2E_PYTHON:-python3}"
T="$(mktemp -d)"

# --- [1/3] generator determinism ---
echo "[1/3] generator determinism"
"$PY" "$HERE/make_synthetic_cohort.py" "$T/g1" >/dev/null
"$PY" "$HERE/make_synthetic_cohort.py" "$T/g2" >/dev/null
for f in data/raw/apa_after_combat.tsv data/raw/CHB.all.meta truth.json genotype/geno.dosage.tsv; do
  cmp -s "$T/g1/$f" "$T/g2/$f" || { echo "FAIL: $f differs between runs"; exit 1; }
done
echo "  PASS"

# --- [2/3] end-to-end runnable chain ---
echo "[2/3] end-to-end runnable chain"
APA_E2E_TMP="$T/e2e" bash "$HERE/run_e2e.sh"

# --- [3/3] planted signal ---
echo "[3/3] planted signal recovered (checked inside run_e2e.sh; log: $T/e2e/logs/planted_signal.log)"
grep -q PASS "$T/e2e/logs/planted_signal.log"
echo "  PASS"
echo "ALL PASS"