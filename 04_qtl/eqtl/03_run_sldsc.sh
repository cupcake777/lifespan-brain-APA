#!/usr/bin/env bash
# Joint apaQTL + eQTL maxCPP S-LDSC: sites, annotations, LD scores, GWAS munging, regressions.
# Work bundle: data/04_qtl/sldsc/; summary tables: data/04_qtl/; figures: figures/04_qtl/.
#SBATCH -J apa_eqtl_sldsc_03
#SBATCH -p cu
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"

# --- config ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
D="$APA_ROOT/data/04_qtl"
W="$D/sldsc"
APA_DIR="${APA_DIR:-$D}"
EQTL_DIR="${EQTL_DIR:-$D/eqtl}"
SLDSC_REF="${SLDSC_REF:-$APA_ROOT/ref/sldsc}"
GWAS_DIR="${GWAS_DIR:-$APA_ROOT/ref/gwas}"
LDSC_DIR="${LDSC_DIR:-$APA_ROOT/ref/ldsc}"
LDSC_PY="${LDSC_PY:-$APA_HPC_HOME/.local/share/mamba/envs/ldsc/bin/python}"
BASE_PY="${BASE_PY:-$APA_HPC_HOME/.local/share/mamba/bin/python3}"
RSCRIPT="${RSCRIPT:-Rscript}"
MERGE_ALLELES_DIR="${MERGE_ALLELES_DIR:-$W/merge_alleles}"
RUN_LOCAL="${RUN_LOCAL:-0}"
RUN_PLOTS="${RUN_PLOTS:-1}"
CHROMOSOMES="${CHROMOSOMES:-$(seq 1 22)}"
STAGES="${STAGES:-prenatal postnatal}"
TRAITS="${TRAITS:-EAS_AD EAS_AIS EAS_ALS EAS_BIP EAS_EP EAS_PD EAS_SCZ EUR_AD EUR_ADHD EUR_ANX EUR_ASD EUR_BIP EUR_EP EUR_INS EUR_INT EUR_MDD EUR_PD EUR_PTSD EUR_SCZ}"
export SCRIPT_DIR W D APA_DIR EQTL_DIR SLDSC_REF GWAS_DIR LDSC_DIR LDSC_PY BASE_PY RSCRIPT MERGE_ALLELES_DIR CHROMOSOMES STAGES TRAITS RUN_LOCAL RUN_PLOTS
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

if [[ -z "${SLURM_JOB_ID:-}" && "${RUN_LOCAL}" != "1" ]]; then
  mkdir -p "$APA_ROOT/logs"
  sbatch -o "$APA_ROOT/logs/03_run_sldsc_%j.out" -e "$APA_ROOT/logs/03_run_sldsc_%j.err" "$0"
  exit 0
fi
mkdir -p "$W/qtl_sites" "$W/ldscores/EAS" "$W/ldscores/EUR" "$W/munged_gwas" "$W/logs" "$D/annot_sd_parts"

# --- maxCPP sites from fine-mapping ---
for stage in ${STAGES}; do
  "$BASE_PY" "${SCRIPT_DIR}/03-1_prep.py" extract "${APA_DIR}/susie/${stage}.SuSiE_summary.parquet" "$W/qtl_sites/${stage}.apaQTL.maxcpp.tsv"
  "$BASE_PY" "${SCRIPT_DIR}/03-1_prep.py" extract "${EQTL_DIR}/susie/${stage}.SuSiE_summary.parquet" "$W/qtl_sites/${stage}.eQTL.maxcpp.tsv"
done

# --- ancestry-matched annotations and ld scores ---
for pop in EAS EUR; do
  for stage in ${STAGES}; do
    for chr in ${CHROMOSOMES}; do
      printf '%s %s %s\n' "$pop" "$stage" "$chr"
    done
  done
done | xargs -r -n 3 -P 4 bash -c '
  set -euo pipefail
  pop=$1; stage=$2; chr=$3
  out="${W}/ldscores/${pop}/${stage}_maxCPP.${chr}"
  "$BASE_PY" "${SCRIPT_DIR}/03-1_prep.py" annot "${SLDSC_REF}/${pop}_hg38/plink_files/1000G.${pop}.hg38.${chr}.bim" \
    "${W}/qtl_sites/${stage}.apaQTL.maxcpp.tsv" "${W}/qtl_sites/${stage}.eQTL.maxcpp.tsv" "${out}.annot.gz"
  "$LDSC_PY" "${LDSC_DIR}/ldsc.py" --l2 \
    --bfile "${SLDSC_REF}/${pop}_hg38/plink_files/1000G.${pop}.hg38.${chr}" \
    --ld-wind-cm 1 --annot "${out}.annot.gz" --print-snps "${SLDSC_REF}/hm3_no_MHC.list.txt" \
    --out "$out" > "${W}/logs/l2.${pop}.${stage}_maxCPP.${chr}.log" 2>&1
' _

# --- common-snp count and annotation sd ---
for pop in EAS EUR; do
  for stage in ${STAGES}; do
    "$BASE_PY" "${SCRIPT_DIR}/03-2_calc_annot.py" "$pop" "$stage" "$SLDSC_REF" "$W/ldscores" "$D/annot_sd_parts/${pop}_${stage}.tsv"
  done
done
awk 'FNR == 1 && NR != 1 {next} {print}' "$D"/annot_sd_parts/*.tsv > "$D/annot_sd.tsv"

# --- gwas allele and sample-size qc ---
for trait in ${TRAITS}; do
  printf '%s\n' "$trait"
done | xargs -r -n 1 -P 4 bash -c '
  set -euo pipefail
  trait=$1; pop=${trait%%_*}
  input="${GWAS_DIR}/${trait}_b19to38_harmonized.tsv"
  if [[ "$trait" == EAS_AD ]]; then input="${GWAS_DIR}/${trait}_b38to38_harmonized.tsv"; fi
  "$LDSC_PY" "${LDSC_DIR}/munge_sumstats.py" --sumstats "$input" --snp SNPID --N-col N \
    --a1 EA --a2 NEA --p P --signed-sumstats BETA,0 \
    --merge-alleles "${MERGE_ALLELES_DIR}/${pop}.hm3_merge_alleles.tsv" \
    --out "${W}/munged_gwas/${trait}" > "${W}/logs/munge.${trait}.log" 2>&1
' _

# --- joint apa-qtl and eqtl regressions ---
for stage in ${STAGES}; do
  for trait in ${TRAITS}; do
    printf '%s %s\n' "$stage" "$trait"
  done
done | xargs -r -n 2 -P 4 bash -c '
  set -euo pipefail
  stage=$1; trait=$2; pop=${trait%%_*}
  out="${W}/results/maxcpp/${pop}/${stage}/${trait}"
  mkdir -p "${out%/*}"
  "$LDSC_PY" "${LDSC_DIR}/ldsc.py" --h2 "${W}/munged_gwas/${trait}.sumstats.gz" \
    --ref-ld-chr "${SLDSC_REF}/${pop}_hg38/baselineLD_v2.2/baselineLD.,${W}/ldscores/${pop}/${stage}_maxCPP." \
    --w-ld-chr "${SLDSC_REF}/${pop}_hg38/weights/weights.hm3_noMHC." \
    --frqfile-chr "${SLDSC_REF}/${pop}_hg38/plink_files/1000G.${pop}.hg38." \
    --overlap-annot --print-coefficients --out "$out" \
    > "${W}/logs/h2.maxcpp.${pop}.${stage}.${trait}.log" 2>&1
' _

# --- manuscript summaries and figures ---
if [[ "${RUN_PLOTS}" == "1" ]]; then
  "${RSCRIPT}" --vanilla "${SCRIPT_DIR}/03-3_plot.R"
fi
echo "MAXCPP_RERUN_COMPLETE"
