#!/usr/bin/env bash
set -euo pipefail

echo "============================================"
echo "refine_bacs_job.sh"
echo "Timestamp : $(date --iso-8601=seconds)"
echo "Host      : $(hostname)"
echo "PWD       : $(pwd)"
echo "RESULTS_DIR   : ${RESULTS_DIR:-}"
echo "REFINE_ENV    : ${REFINE_BACS_ENV:-}"
echo "REFINE_RUN_ARGS: ${REFINE_RUN_ARGS:-}"
echo "CPUS         : ${CPUS:-}"
echo "============================================"

bash workflow/run_refine_annot_smorf_bacs.sh \
  --run ${REFINE_RUN_ARGS} \
  --cpus "${CPUS:-8}"
