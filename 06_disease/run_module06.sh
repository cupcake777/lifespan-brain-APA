#!/usr/bin/env bash
# Module 06 sLDSC: APA-only sLDSC, then joint apaQTL+eQTL maxCPP sLDSC.
# Submits itself to Slurm unless SLURM_JOB_ID is set or RUN_LOCAL=1.
#SBATCH --job-name=module06_sldsc
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=cu

set -euo pipefail

# --- config ---

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
APA_OUT="${APA_ROOT}/data/06_disease"
JOINT_OUT="${APA_ROOT}/data/04_qtl/sldsc"
APA_SCRIPT="${APA_SCRIPT:-${MODULE_DIR}/sldsc/run_sldsc.sbatch}"
JOINT_SCRIPT="${JOINT_SCRIPT:-${MODULE_DIR}/../04_qtl/eqtl/03_run_sldsc.sh}"
RUN_LOCAL="${RUN_LOCAL:-0}"

# --- submit to slurm ---

if [[ -z "${SLURM_JOB_ID:-}" && "${RUN_LOCAL}" != "1" ]]; then
    mkdir -p "${APA_ROOT}/logs"
    export APA_SCRIPT JOINT_SCRIPT RUN_LOCAL
    sbatch --chdir="${MODULE_DIR}" --output="${APA_ROOT}/logs/module06_%j.out" "${BASH_SOURCE[0]}"
    exit 0
fi

mkdir -p "${APA_OUT}" "${JOINT_OUT}"

# --- APA-only sLDSC ---

echo "[module06] APA-only sLDSC"
env \
    SCRIPT_DIR="$(cd "$(dirname "${APA_SCRIPT}")" && pwd)" \
    RUN_PLOTS="${APA_RUN_PLOTS:-${RUN_PLOTS:-1}}" \
    bash "${APA_SCRIPT}"

# --- joint apaQTL+eQTL maxCPP sLDSC ---

echo "[module06] joint apaQTL+eQTL maxCPP sLDSC"
joint_env=(RUN_LOCAL=1)
for name in APA_DIR EQTL_DIR SLDSC_REF GWAS_DIR LDSC_DIR LDSC_PY BASE_PY RSCRIPT MERGE_ALLELES_DIR RUN_PLOTS CHROMOSOMES STAGES TRAITS; do
    var_name="JOINT_${name}"
    value="${!var_name:-}"
    if [[ -n "${value}" ]]; then
        joint_env+=("${name}=${value}")
    fi
done
env "${joint_env[@]}" bash "${JOINT_SCRIPT}"

# --- completion gate ---

apa_result_count="$(find "${APA_OUT}/regression" -type f -name '*.results' | wc -l)"
joint_result_count="$(find "${JOINT_OUT}/results/maxcpp" -type f -name '*.results' | wc -l)"
[[ -s "${APA_OUT}/sldsc_enrichment_summary.tsv" ]]
(( apa_result_count > 0 ))
(( joint_result_count > 0 ))

# --- manifest ---

manifest="${APA_OUT}/module06_manifest.tsv"
printf '%s\n' 'analysis	output_root	estimand' > "${manifest}"
printf '%s\n' "apa_only	${APA_OUT}	APA-only custom annotations" >> "${manifest}"
printf '%s\n' "joint_apa_eqtl_maxcpp	${JOINT_OUT}	stage-specific joint apaQTL/eQTL maxCPP annotations" >> "${manifest}"
echo "[module06] wrote ${manifest}"
