#!/usr/bin/env bash
set -euo pipefail

echo "============================================"
echo "refine_euks_job.sh"
echo "Timestamp : $(date --iso-8601=seconds)"
echo "Host      : $(hostname)"
echo "PWD       : $(pwd)"
echo "RESULTS_DIR         : ${RESULTS_DIR:-}"
echo "REFINE_EUKS_ENV     : ${REFINE_EUKS_ENV:-}"
echo "REFINE_EUKS_RUN_ARGS: ${REFINE_EUKS_RUN_ARGS:-}"
echo "CPUS                : ${CPUS:-}"
echo "SMORFS_WORK_ROOT    : ${SMORFS_WORK_ROOT:-}"
echo "============================================"

bash workflow/run_refine_annot_smorf_euks.sh \
  --refine-env "${REFINE_EUKS_ENV:-refine_annot_smorf_euks_env}" \
  --run ${REFINE_EUKS_RUN_ARGS} \
  --cpus "${CPUS:-8}"
