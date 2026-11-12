#!/usr/bin/env bash
# Checks the publication tree: required files, no symlinks, no stale dated-folder calls, shell syntax.
set -euo pipefail

package_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# --- symlinks ---
if find "${package_root}" -type l -print -quit | grep -q .; then
  echo "ERROR: publication tree contains filesystem symlinks" >&2
  find "${package_root}" -type l -print >&2
  exit 1
fi

# --- required files ---
required=(
  MANIFEST.tsv
  README.md
  00_geno/run_module00.sh
  01_quant/rna_quant.sh
  02_traj/run.sh
  02_traj/README_data.md
  02_traj/README_outputs.md
  03_omics/run_protein.sh
  04_qtl/apaqtl/qtl_pipe.sh
  05_gxs/01_tensorqtl.sh
  05_gxs/plot_fig.R
  06_disease/run_module06.sh
  07_coloc/easycoloc/LICENSE
  07_coloc/easycoloc/environment.yml
  07_coloc/coloc/input.sh
  07_coloc/coloc/06_plot_summary.R
  07_coloc/func/run_analysis.sh
  07_coloc/func/config.yml
)
for path in "${required[@]}"; do
  test -f "${package_root}/${path}" || {
    echo "ERROR: missing ${path}" >&2
    exit 1
  }
done

# --- stale dated-folder calls ---
if grep -R -n -E \
  --include='*.sh' --include='*.sbatch' --include='*.slurm' \
  --include='*.R' --include='*.r' --include='*.py' \
  --exclude='check_portable_tree.sh' \
  'pub/1010|\.\./\.\./1010' "${package_root}"; then
  echo "ERROR: publication scripts still call the retained 1010 tree" >&2
  exit 1
fi

# --- shell syntax ---
while IFS= read -r script; do
  bash -n "${script}"
done < <(find "${package_root}" -type f \( -name '*.sh' -o -name '*.sbatch' -o -name '*.slurm' \) -print)

# --- manifest columns ---
awk -F '\t' 'NF != 13 {exit 1}' "${package_root}/MANIFEST.tsv"
echo "Portable publication tree check passed."
