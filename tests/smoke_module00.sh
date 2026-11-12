#!/usr/bin/env bash
# Smoke test for module 00: joint-call merge and stage-specific variant filter.
# Stubs tabix, bcftools and plink2; all outputs go to a scratch APA_ROOT.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
module="${root}/00_geno"
tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT
export APA_ROOT="${tmp}/root"
export TMPDIR="${tmp}/work"

stub_bin="${tmp}/bin"
mkdir -p "${stub_bin}"

cat > "${stub_bin}/tabix" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
touch "${@: -1}.tbi"
EOF

cat > "${stub_bin}/bcftools" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
if [[ "$1" != concat ]]; then
  exit 0
fi
out=""
list=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) out="$2"; shift 2 ;;
    -f) list="$2"; shift 2 ;;
    *) shift ;;
  esac
done
[[ -n "$out" && -n "$list" ]]
cp "$list" "$out"
EOF

cat > "${stub_bin}/plink2" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
out=""
mode=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --out) out="$2"; shift 2 ;;
    --write-snplist) mode=snplist; shift ;;
    --make-pgen) mode=pgen; shift ;;
    *) shift ;;
  esac
done
[[ -n "$out" ]]
case "$mode" in
  snplist)
    if [[ "$out" == *pre_qc ]]; then
      printf '%s\n' '1:100:A:G' '1:200:C:T' > "${out}.snplist"
    else
      printf '%s\n' '1:200:C:T' '1:300:G:A' > "${out}.snplist"
    fi
    ;;
  pgen)
    : > "${out}.pgen"
    : > "${out}.pvar"
    : > "${out}.psam"
    ;;
  *) exit 1 ;;
esac
EOF
chmod +x "${stub_bin}"/*

# --- merge ordering ---
merge_dir="${APA_ROOT}/data/00_geno/joint_calling"
mkdir -p "${merge_dir}"
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY; do
  printf '%s\n' "${chr}" > "${merge_dir}/joint_genotyped.${chr}.vcf.gz"
done
PATH="${stub_bin}:${PATH}" bash "${module}/02_merge_clean.sh" >/dev/null
expected="${tmp}/expected.list"
printf '%s\n' chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
  | sed 's#^#joint_genotyped.#; s#$#.vcf.gz#' > "${expected}"
diff -u "${expected}" "${merge_dir}/vcf_list.txt"
[[ -e "${APA_ROOT}/data/00_geno/joint_all_chr.vcf.gz.tbi" ]]

# --- stage intersection ---
meta="${tmp}/meta.tsv"
cat > "${meta}" <<'EOF'
sample	batch	RIN	sex	age	LifeStage
pre1	b1	7	Female	1	Prenatal
post1	b1	7	Male	10	Postnatal
EOF
vcf="${tmp}/input.vcf.gz"
: > "${vcf}"
PATH="${stub_bin}:${PATH}" RAW_VCF="${vcf}" META_FILE="${meta}" bash "${module}/00_4_stage_filter.sh" >/dev/null
stage_dir=("${TMPDIR}"/stage_qc_*)
printf '%s\n' 'chr1:200:C:T' > "${tmp}/expected.common"
diff -u "${tmp}/expected.common" "${stage_dir[0]}/common_variants.txt"
[[ -e "${APA_ROOT}/data/00_geno/genotype_clean.pgen" ]]

printf '%s\n' 'Module 00 smoke test passed.'
