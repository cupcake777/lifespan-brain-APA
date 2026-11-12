#!/usr/bin/env bash
# Synthetic recovery test for the S2 element tables (03_omics/s2_element.py). Exit 0 on success.
# Inputs are generated into a temporary folder; nothing is written under release/.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${PYTHON:-python3}"
fail=0

# --- run ---
echo "== S2 element tables: synthetic recovery"
if "${py}" "${here}/test_s2_recovery.py"; then
  echo "PASS s2_recovery"
else
  echo "FAIL s2_recovery"; fail=1
fi

# --- report ---
if [[ "${fail}" -eq 0 ]]; then
  echo "ALL OMICS SYNTHETIC TESTS PASSED"
else
  echo "OMICS SYNTHETIC TESTS FAILED"
fi
exit "${fail}"
