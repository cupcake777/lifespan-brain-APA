#!/usr/bin/env bash
# Run the runnable part of the module 02 chain (01_cov to 04_cluster) on a generated synthetic cohort.
# Steps that need unavailable tools are reported as skipped; exit non-zero only if a runnable step fails.
# Env: APA_E2E_TMP (work dir), APA_E2E_RSCRIPT (Rscript for 04), APA_E2E_PYTHON (python with numpy, pandas).
set -uo pipefail

# --- setup ---
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REL="$(cd "$HERE/../../.." && pwd)"
TMP="${APA_E2E_TMP:-$(mktemp -d)}"
COH="$TMP/cohort"; RUN="$TMP/run"; LOG="$TMP/logs"
mkdir -p "$COH" "$RUN/pipeline" "$RUN/data/raw" "$RUN/outputs/intermediates" "$LOG"
STATUS="$TMP/step_status.tsv"; printf 'step\tscript\tstatus\tdetail\n' > "$STATUS"
fail=0
PY="${APA_E2E_PYTHON:-python3}"

# --- synthetic cohort and module 02 scripts ---
"$PY" "$HERE/make_synthetic_cohort.py" "$COH" > "$LOG/cohort.log" 2>&1 || { echo "cohort generation failed"; exit 1; }
cp "$COH/data/raw/"* "$RUN/data/raw/"
cp "$REL/02_traj/"[0-9]*.R "$RUN/pipeline/"
cp "$COH/truth.json" "$TMP/truth.json"

# --- runnable steps ---
run_step () {  # name script [env...]
  local name="$1" script="$2"; shift 2
  ( cd "$RUN" && env "$@" timeout 1500 Rscript "pipeline/$script" ) > "$LOG/$name.log" 2>&1
  local rc=$?
  if [ $rc -eq 0 ]; then printf '%s\t%s\tran\t\n' "$name" "$script" >> "$STATUS"
  else
    local err; err="$(grep -m1 -E 'Error|error' "$LOG/$name.log" | tr '\t' ' ')"
    printf '%s\t%s\tfailed\trc=%s %s\n' "$name" "$script" "$rc" "$err" >> "$STATUS"; fail=1
  fi
}

run_step 01_cov 01_cov.R
run_step 02_support 02_support.R
run_step 02_gam 02_gam.R
run_step 02_prediction 02_prediction.R
run_step 03_dynamic 03_dynamic.R APA_DYNAMIC_N_PERM=50 APA_DYNAMIC_SENS_N_PERM=20 APA_DYNAMIC_SITES_PER=60

# --- 04_cluster (needs R package bluster) ---
CLUSTER_R="${APA_E2E_RSCRIPT:-Rscript}"
if "$CLUSTER_R" -e 'quit(status = if (requireNamespace("bluster", quietly = TRUE)) 0 else 1)' >/dev/null 2>&1; then
  ( cd "$RUN" && APA_CLUSTER_BOOT_ITER=20 timeout 1500 "$CLUSTER_R" pipeline/04_cluster.R ) > "$LOG/04_cluster.log" 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then printf '04_cluster\t04_cluster.R\tran\t\n' >> "$STATUS"
  else printf '04_cluster\t04_cluster.R\tfailed\trc=%s\n' "$rc" >> "$STATUS"; fail=1; fi
else
  printf '04_cluster\t04_cluster.R\tskipped\tneeds R package bluster (Bioconductor); not in the default r environment\n' >> "$STATUS"
fi

# --- planted-signal check ---
if [ $fail -eq 0 ]; then
  "$PY" "$HERE/check_planted_signal.py" "$RUN" "$TMP/truth.json" > "$LOG/planted_signal.log" 2>&1 \
    && printf 'planted_signal\tcheck_planted_signal.py\tran\t\n' >> "$STATUS" \
    || { printf 'planted_signal\tcheck_planted_signal.py\tfailed\tsee logs\n' >> "$STATUS"; fail=1; }
fi

cat "$STATUS"
echo "work dir: $TMP"
exit $fail