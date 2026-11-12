#!/usr/bin/env bash
# Smoke test for module 06: APA-only sLDSC and joint apaQTL+eQTL maxCPP sLDSC, with stub tools.
# All outputs go to a scratch APA_ROOT.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT
export APA_ROOT="${tmp}/root"
export APA_HPC_HOME="${tmp}/hpc"
export APA_GROUP_SHARE="${tmp}/group_share"

# --- stub tools ---
stub_bin="${tmp}/bin"
mkdir -p "${stub_bin}"

cat > "${stub_bin}/mock_python" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
script="$(basename "$1")"
shift
case "${script}" in
  03_annot.py)
    mkdir -p "$(dirname "$4")"
    printf '%s\n' 'CHR BP SNP CM annotation' > "$4"
    ;;
  05_format_gwas.py)
    mkdir -p "$(dirname "$2")"
    printf '%s\n' 'SNP A1 A2 N Z FRQ' > "$2"
    ;;
  03-1_prep.py)
    mode="$1"
    if [[ "${mode}" == extract ]]; then
      out="$3"
      mkdir -p "$(dirname "${out}")"
      printf '%s\n' '1\t100\tA\tG\t0.9' > "${out}"
    elif [[ "${mode}" == annot ]]; then
      out="$5"
      mkdir -p "$(dirname "${out}")"
      printf '%s\n' 'CHR\tBP\tSNP\tCM\tapaQTL\teQTL' > "${out}"
    else
      exit 1
    fi
    ;;
  03-2_calc_annot.py)
    out="$5"
    mkdir -p "$(dirname "${out}")"
    printf '%s\n' 'population\tstage\tannot_set\tannotation\tn_snps\tmean\tsd' > "${out}"
    printf '%s\n' 'mock\tmock\tmaxCPP\tapaQTL\t1\t0.1\t0.1' >> "${out}"
    printf '%s\n' 'mock\tmock\tmaxCPP\teQTL\t1\t0.1\t0.1' >> "${out}"
    ;;
  *)
    exit 1
    ;;
esac
EOF

cat > "${stub_bin}/mock_ldsc_python" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
script="$(basename "$1")"
shift
if [[ "${script}" == "07_collect.py" ]]; then
  output="$2"
  mkdir -p "$(dirname "${output}")"
  printf '%s\n' 'Trait\tPopulation\tAnnotation' > "${output}"
  exit 0
fi
out=""
is_l2=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --l2) is_l2=1; shift ;;
    --out) out="$2"; shift 2 ;;
    *) shift ;;
  esac
done

case "${script}" in
  ldsc.py)
    [[ -n "${out}" ]]
    mkdir -p "$(dirname "${out}")"
    if [[ "${is_l2}" == "1" ]]; then
      printf '%s\n' 'L2 mock' > "${out}.l2.ldscore.gz"
    else
      printf '%s\n' 'Category\tCoefficient\tCoefficient_std_error\tCoefficient_z-score\tEnrichment\tEnrichment_std_error\tEnrichment_p' \
        'mock\t0.1\t0.01\t10\t1.1\t0.1\t0.001' > "${out}.results"
    fi
    ;;
  munge_sumstats.py)
    [[ -n "${out}" ]]
    mkdir -p "$(dirname "${out}")"
    printf '%s\n' 'SNP\tA1\tA2\tN\tZ\tFRQ' > "${out}.sumstats.gz"
    ;;
  08_plot_sldsc.py)
    ;;
  *)
    exit 1
    ;;
esac
EOF

chmod +x "${stub_bin}/mock_python" "${stub_bin}/mock_ldsc_python"

# --- scratch dirs ---
run_root="${tmp}/run"
harm_dir="${tmp}/harmonized"
ref_root="${tmp}/ref"
joint_gwas="${tmp}/joint_gwas"
merge_dir="${tmp}/merge_alleles"
mkdir -p "${harm_dir}" "${ref_root}" "${joint_gwas}" "${merge_dir}"

# --- synthetic inputs ---
for annotation in prenatal_apa_qtl postnatal_apa_qtl interaction_apa_qtl lifespan_apa_qtl; do
  mkdir -p "${APA_ROOT}/data/06_disease/qtl_annot"
  printf '%s\n' 'chr1 99 101' > "${APA_ROOT}/data/06_disease/qtl_annot/${annotation}.bed"
done
printf '%s\n' 'SNPID\tvariant_id\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN' > "${harm_dir}/EUR_MDD_b19to38_harmonized.tsv"
printf '%s\n' 'SNPID\tvariant_id\tCHR\tPOS\tEA\tNEA\tEAF\tBETA\tSE\tP\tN' > "${harm_dir}/EAS_AD_b38to38_harmonized.tsv"
printf '%s\n' '1 1 rs1 0 100 A G' > "${ref_root}/EUR.1.bim"
printf '%s\n' '1 1 rs1 0 100 A G' > "${ref_root}/EAS.1.bim"
for pop in EAS EUR; do
  printf '%s\n' 'SNP' > "${merge_dir}/${pop}.hm3_merge_alleles.tsv"
  printf '%s\n' 'SNP MAF' > "${ref_root}/${pop}.frq"
done

for stage in prenatal postnatal; do
  mkdir -p "${tmp}/apa_input/susie" "${tmp}/eqtl_input/susie"
  : > "${tmp}/apa_input/susie/${stage}.SuSiE_summary.parquet"
  : > "${tmp}/eqtl_input/susie/${stage}.SuSiE_summary.parquet"
done
for trait in EAS_AD EUR_MDD; do
  : > "${joint_gwas}/${trait}_b19to38_harmonized.tsv"
done
: > "${joint_gwas}/EAS_AD_b38to38_harmonized.tsv"

# --- run module 06 ---
RUN_LOCAL=1 \
APA_RUN_PLOTS=0 \
JOINT_RUN_PLOTS=0 \
CHROMOSOMES=1 \
EUR_TRAITS_OVERRIDE='MDD' \
EAS_TRAITS_OVERRIDE='AD' \
PYTHON_BIN="${stub_bin}/mock_python" \
LDSC_PYTHON="${stub_bin}/mock_ldsc_python" \
LDSC_BIN="${tmp}/ldsc.py" \
GWAS_HARM_DIR="${harm_dir}" \
EUR_BFILE="${ref_root}/EUR" \
EAS_BFILE="${ref_root}/EAS" \
HM3="${ref_root}/hm3.list" \
EUR_BASELINE="${ref_root}/eur.baseline." \
EAS_BASELINE="${ref_root}/eas.baseline." \
EUR_WEIGHTS="${ref_root}/eur.weights." \
EAS_WEIGHTS="${ref_root}/eas.weights." \
EUR_FRQFILE="${ref_root}/eur.frq." \
EAS_FRQFILE="${ref_root}/eas.frq." \
JOINT_APA_DIR="${tmp}/apa_input" \
JOINT_EQTL_DIR="${tmp}/eqtl_input" \
JOINT_SLDSC_REF="${ref_root}" \
JOINT_GWAS_DIR="${joint_gwas}" \
JOINT_LDSC_DIR="${tmp}/ldsc" \
JOINT_LDSC_PY="${stub_bin}/mock_ldsc_python" \
JOINT_BASE_PY="${stub_bin}/mock_python" \
JOINT_MERGE_ALLELES_DIR="${merge_dir}" \
JOINT_CHROMOSOMES=1 \
JOINT_STAGES='prenatal postnatal' \
JOINT_TRAITS='EAS_AD EUR_MDD' \
PATH="${stub_bin}:${PATH}" \
bash "${root}/06_disease/run_module06.sh" > "${tmp}/module06.log"

# --- checks ---
apa_out="${APA_ROOT}/data/06_disease"
joint_out="${APA_ROOT}/data/04_qtl/sldsc"
[[ -s "${apa_out}/module06_manifest.tsv" ]]
[[ "$(find "${apa_out}/regression" -name '*.results' | wc -l)" -eq 8 ]]
[[ -s "${apa_out}/sldsc_enrichment_summary.tsv" ]]
[[ "$(find "${joint_out}/results/maxcpp" -name '*.results' | wc -l)" -eq 4 ]]
[[ "$(find "${joint_out}/ldscores" -name '*.l2.ldscore.gz' | wc -l)" -eq 4 ]]
awk -F '\t' 'NR == 2 && $1 == "apa_only" {apa=1} NR == 3 && $1 == "joint_apa_eqtl_maxcpp" {joint=1} END {exit !(apa && joint)}' "${apa_out}/module06_manifest.tsv"
grep -q -E '^\[module06\] APA-only sLDSC$' "${tmp}/module06.log"
grep -q -E '^\[module06\] joint apaQTL\+eQTL maxCPP sLDSC$' "${tmp}/module06.log"

printf '%s\n' 'Module 06 smoke test passed.'
