#!/usr/bin/env bash
# Synthetic tests for the disease track (modules 06 and 07); no real data are read.
# PYTHON needs numpy and pandas; RSCRIPT needs coloc 5.x and susieR.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"
export RSCRIPT="${RSCRIPT:-Rscript}"
export PYTHONDONTWRITEBYTECODE=1
"$PY" "$HERE/test_sldsc_2stage.py"
"$PY" "$HERE/test_magma.py"
"$PY" "$HERE/test_coloc.py"
echo "ALL DISEASE SYNTHETIC TESTS PASSED"
