#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${repo_root}"

bash -n easycoloc
bash -n tools/run_pipeline_managed.sh
bash -n tools/watch_output_dir.sh

tmp_dir="$(mktemp -d)"
trap 'rm -rf "${tmp_dir}"' EXIT

./easycoloc init "${tmp_dir}/project"

test -f "${tmp_dir}/project/config/global.yml"
test -f "${tmp_dir}/project/config/gwas.yml"
test -f "${tmp_dir}/project/config/qtl.yml"

./easycoloc refs \
  --global "${tmp_dir}/project/config/global.yml" \
  --gwas "${tmp_dir}/project/config/gwas.yml" \
  --qtl "${tmp_dir}/project/config/qtl.yml" >/dev/null

./easycoloc help | grep "report-web" >/dev/null
./easycoloc help | grep "harmonize" >/dev/null
./easycoloc help | grep "quick" >/dev/null
./easycoloc quick --help | grep -- "--gwas-id ID --qtl-id ID" >/dev/null
./easycoloc quick --list | grep "example_eqtl" >/dev/null
./easycoloc harmonize --help >/dev/null
./easycoloc report-web --help >/dev/null
./easycoloc run --emit-slurm-array-template --shard-count 3 --output-dir "${tmp_dir}/project/results" | grep "#SBATCH --array=1-3" >/dev/null

echo "[SMOKE] cli smoke test passed"
