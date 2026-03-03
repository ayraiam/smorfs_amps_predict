#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------
# Use the SAME env creation logic as refine stage
# (single source of truth)
# --------------------------------------------------

ENV_NAME="${REFINE_BACS_ENV:-refine_annot_smorf_bacs_env}"
ENV_PREFIX="envs/${ENV_NAME}"

# Ensure env exists using the official helper (idempotent)
bash workflow/run_refine_annot_smorf_bacs.sh --create-env --refine-env "${ENV_NAME}"

# Activate it
conda activate "${ENV_PREFIX}"

# Sanity check
command -v mmseqs >/dev/null 2>&1 || {
  echo "ERROR: mmseqs not found after activating ${ENV_PREFIX}" >&2
  exit 1
}

echo "[INFO] Using mmseqs from: $(which mmseqs)"

echo "============================================"
echo "mmseqs_global_job.sh"
echo "Timestamp : $(date --iso-8601=seconds)"
echo "Host      : $(hostname)"
echo "PWD       : $(pwd)"
echo "RESULTS_DIR : ${RESULTS_DIR:-}"
echo "THREADS    : ${THREADS:-}"
echo "MIN_SEQ_ID : ${MIN_SEQ_ID:-}"
echo "COV        : ${COV:-}"
echo "COV_MODE   : ${COV_MODE:-}"
echo "============================================"

bash workflow/run_mmseqs_global_cluster.sh \
  --results-dir "${RESULTS_DIR:-results}" \
  --threads "${THREADS:-8}" \
  --min-seq-id "${MIN_SEQ_ID:-0.95}" \
  --cov "${COV:-0.8}" \
  --cov-mode "${COV_MODE:-1}"
