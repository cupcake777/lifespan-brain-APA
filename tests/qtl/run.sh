#!/usr/bin/env bash
# Synthetic QTL tests: simulated inputs only, nothing written under release/.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RELEASE="$(cd "$HERE/../.." && pwd)"
PY="${PYTHON:-python3}"
RS="${RSCRIPT:-Rscript}"

# --- run ---
"$PY" "$HERE/test_hc3_refit.py"
"$PY" "$HERE/test_gating.py"
"$RS" "$HERE/test_mediation.R" "$RELEASE/04_qtl/eqtl/04-1_mediate.R"
echo "ALL QTL SYNTHETIC TESTS PASSED"
