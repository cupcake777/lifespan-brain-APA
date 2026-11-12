#!/usr/bin/env bash
# Submit full PDUI TWAS model training as a SLURM array, then a combine job.
# Logs go to $APA_ROOT/logs/training; the model tree stays under ROOT.

# --- config ---
STAGE="${STAGE:-prenatal}"
CHUNK_SIZE="${CHUNK_SIZE:-200}"
THREADS="${THREADS:-1}"
N_FOLDS="${N_FOLDS:-5}"
MEM="${MEM:-8G}"
TIME="${TIME:-48:00:00}"
ARRAY_CONCURRENCY="${ARRAY_CONCURRENCY:-24}"
COMBINE_MEM="${COMBINE_MEM:-16G}"
COMBINE_TIME="${COMBINE_TIME:-12:00:00}"

set -euo pipefail

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# --- paths ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT="${APA_ROOT}/data/06_disease/twas"
LOG_DIR="${APA_ROOT}/logs/training/${STAGE}_array_$(date +%Y%m%d_%H%M%S)"
mkdir -p "${LOG_DIR}"

# --- manifest and array size ---
python "${SCRIPT_DIR}/02_twas_train.py" manifest \
  --stage "${STAGE}" \
  --chunk-size "${CHUNK_SIZE}" > "${LOG_DIR}/${STAGE}_manifest.tsv"

ARRAY_SIZE="$(awk 'NR > 1 {n++} END {print n + 0}' "${ROOT}/models/${STAGE}/chunk_manifest.tsv")"
if [[ "${ARRAY_SIZE}" -lt 1 ]]; then
  echo "ERROR: no chunks found for ${STAGE}" >&2
  exit 1
fi

# --- submit array and combine ---
ARRAY_JOB_RAW="$(
  sbatch \
    --parsable \
    --job-name="twas_${STAGE}_chunk" \
    --array="1-${ARRAY_SIZE}%${ARRAY_CONCURRENCY}" \
    --cpus-per-task="${THREADS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --output="${LOG_DIR}/${STAGE}_chunk_%A_%a.out" \
    --error="${LOG_DIR}/${STAGE}_chunk_%A_%a.err" \
    --wrap="cd '${ROOT}' && export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 && python '${SCRIPT_DIR}/02_twas_train.py' train-chunk --stage '${STAGE}' --chunk-id \"\${SLURM_ARRAY_TASK_ID}\" --chunk-size '${CHUNK_SIZE}' --threads '${THREADS}' --n-folds '${N_FOLDS}' --clean-tmp"
)"
ARRAY_JOB_ID="${ARRAY_JOB_RAW%%;*}"

COMBINE_JOB_RAW="$(
  sbatch \
    --parsable \
    --dependency="afterok:${ARRAY_JOB_ID}" \
    --job-name="twas_${STAGE}_combine" \
    --cpus-per-task=1 \
    --mem="${COMBINE_MEM}" \
    --time="${COMBINE_TIME}" \
    --output="${LOG_DIR}/${STAGE}_combine_%j.out" \
    --error="${LOG_DIR}/${STAGE}_combine_%j.err" \
    --wrap="cd '${ROOT}' && python '${SCRIPT_DIR}/02_twas_train.py' combine --stage '${STAGE}'"
)"
COMBINE_JOB_ID="${COMBINE_JOB_RAW%%;*}"

cat <<EOF
stage	${STAGE}
array_job_id	${ARRAY_JOB_ID}
combine_job_id	${COMBINE_JOB_ID}
array_size	${ARRAY_SIZE}
array_concurrency	${ARRAY_CONCURRENCY}
chunk_size	${CHUNK_SIZE}
threads_per_task	${THREADS}
log_dir	${LOG_DIR}
manifest	${ROOT}/models/${STAGE}/chunk_manifest.tsv
expected_db	${ROOT}/models/${STAGE}/${STAGE}.full_twas.db
expected_cov	${ROOT}/models/${STAGE}/${STAGE}.full_twas.cov
expected_detail	${ROOT}/models/${STAGE}/${STAGE}.training_detail.tsv
expected_combine_summary	${ROOT}/models/${STAGE}/combine_summary.tsv
EOF
