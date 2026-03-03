#!/usr/bin/env bash
set -euo pipefail

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
