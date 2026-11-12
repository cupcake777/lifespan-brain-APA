#!/usr/bin/env bash
# S-LDSC on the G x stage maxCPP set: baselineLD plus the interaction annotation only.
# Work bundle data/04_qtl/interaction_sldsc/; munged GWAS and merge alleles in data/04_qtl/sldsc/.
#SBATCH -J apa_interaction_sldsc_04
#SBATCH -p cu
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
set -euo pipefail
: "${APA_ROOT:?source config/paths.sh}"
: "${APA_HPC_HOME:?source config/paths.sh}"

SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
D="$APA_ROOT/data/04_qtl"
W="$D/interaction_sldsc"
SLDSC_REF="$APA_HPC_HOME/project/.ref/sLDSC"
GWAS_DIR="$APA_HPC_HOME/project/.ref/GWAS"
LDSC_DIR="$APA_HPC_HOME/project/work/ref_repo/ldsc"
LDSC_PY="$APA_HPC_HOME/.local/share/mamba/envs/ldsc/bin/python"
BASE_PY="$APA_HPC_HOME/.local/share/mamba/bin/python3"
RSCRIPT="$APA_HPC_HOME/.local/share/mamba/bin/Rscript"
TRAITS="EAS_AD EAS_AIS EAS_ALS EAS_BIP EAS_EP EAS_PD EAS_SCZ EUR_AD EUR_ADHD EUR_ANX EUR_ASD EUR_BIP EUR_EP EUR_INS EUR_INT EUR_MDD EUR_PD EUR_PTSD EUR_SCZ"
export SCRIPT_DIR D W SLDSC_REF GWAS_DIR LDSC_DIR LDSC_PY BASE_PY RSCRIPT
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  mkdir -p "$APA_ROOT/logs"
  sbatch -o "$APA_ROOT/logs/04-3_sldsc_interaction_%j.out" -e "$APA_ROOT/logs/04-3_sldsc_interaction_%j.err" "$0"
  exit 0
fi
mkdir -p "$W/qtl_sites" "$W/ldscores/EAS" "$W/ldscores/EUR" "$W/logs" "$W/results/interaction_maxcpp" "$D/annot_sd_parts"

# --- merge per-chr max pip ---
"$BASE_PY" "$SCRIPT_DIR/04-3_prep_inter.py" merge "$W/qtl_sites/interaction.maxcpp.tsv"

# --- ancestry-matched annotations and ld scores ---
for pop in EAS EUR; do
  for chr in {1..22}; do printf '%s %s\n' "$pop" "$chr"; done
done | xargs -r -n 2 -P 4 bash -c '
  set -euo pipefail
  pop=$1; chr=$2
  out="$W/ldscores/${pop}/interaction_maxCPP.${chr}"
  "$BASE_PY" "$SCRIPT_DIR/04-3_prep_inter.py" annot "$SLDSC_REF/${pop}_hg38/plink_files/1000G.${pop}.hg38.${chr}.bim" \
    "$W/qtl_sites/interaction.maxcpp.tsv" "${out}.annot.gz"
  "$LDSC_PY" "$LDSC_DIR/ldsc.py" --l2 --bfile "$SLDSC_REF/${pop}_hg38/plink_files/1000G.${pop}.hg38.${chr}" \
    --ld-wind-cm 1 --annot "${out}.annot.gz" --print-snps "$SLDSC_REF/hm3_no_MHC.list.txt" \
    --out "$out" > "$W/logs/l2.${pop}.interaction_maxCPP.${chr}.log" 2>&1
' _

# --- common-snp count and annotation sd ---
for pop in EAS EUR; do
  "$BASE_PY" "$SCRIPT_DIR/04-3_prep_inter.py" calcsd "$pop" "$SLDSC_REF" "$W/ldscores" "$D/annot_sd_parts/${pop}_interaction.tsv"
done
awk 'FNR == 1 && NR != 1 {next} {print}' "$D/annot_sd_parts/EAS_interaction.tsv" "$D/annot_sd_parts/EUR_interaction.tsv" > "$D/annot_sd_interaction.tsv"

# --- reuse munged gwas; munge only if missing ---
for trait in $TRAITS; do
  if [[ ! -s "$D/sldsc/munged_gwas/${trait}.sumstats.gz" ]]; then
    pop=${trait%%_*}
    input="${GWAS_DIR}/${trait}_b19to38_harmonized.tsv"
    if [[ "$trait" == EAS_AD ]]; then input="${GWAS_DIR}/${trait}_b38to38_harmonized.tsv"; fi
    "$LDSC_PY" "${LDSC_DIR}/munge_sumstats.py" --sumstats "$input" --snp SNPID --N-col N \
      --a1 EA --a2 NEA --p P --signed-sumstats BETA,0 \
      --merge-alleles "$D/sldsc/merge_alleles/${pop}.hm3_merge_alleles.tsv" \
      --out "$D/sldsc/munged_gwas/${trait}" > "$W/logs/munge.${trait}.log" 2>&1
  fi
done

# --- interaction regressions: baselineLD + interaction_maxCPP only ---
for trait in $TRAITS; do printf '%s\n' "$trait"; done | xargs -r -n 1 -P 4 bash -c '
  set -euo pipefail
  trait=$1; pop=${trait%%_*}
  out="$W/results/interaction_maxcpp/${pop}/${trait}"
  mkdir -p "${out%/*}"
  "$LDSC_PY" "$LDSC_DIR/ldsc.py" --h2 "$D/sldsc/munged_gwas/${trait}.sumstats.gz" \
    --ref-ld-chr "$SLDSC_REF/${pop}_hg38/baselineLD_v2.2/baselineLD.,$W/ldscores/${pop}/interaction_maxCPP." \
    --w-ld-chr "$SLDSC_REF/${pop}_hg38/weights/weights.hm3_noMHC." \
    --frqfile-chr "$SLDSC_REF/${pop}_hg38/plink_files/1000G.${pop}.hg38." \
    --overlap-annot --print-coefficients --out "$out" > "$W/logs/h2.interaction_maxcpp.${pop}.${trait}.log" 2>&1
' _

# --- collect tau star summary ---
"$RSCRIPT" --vanilla "$SCRIPT_DIR/04-4_collect_inter.R"
echo "INTERACTION_SLDSC_COMPLETE"
