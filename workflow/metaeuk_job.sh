#!/usr/bin/env bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

echo "============================================"
echo "metaeuk_job.sh"
echo "Timestamp        : $(date --iso-8601=seconds)"
echo "Host             : $(hostname)"
echo "PWD              : $(pwd)"
echo "RESULTS_DIR      : ${RESULTS_DIR:-}"
echo "SMORFS_WORK_ROOT : ${SMORFS_WORK_ROOT:-}"
echo "METAEUK_DB       : ${METAEUK_DB:-}"
echo "METAEUK_RUN_ARGS : ${METAEUK_RUN_ARGS:-}"
echo "CPUS             : ${CPUS:-}"
echo "============================================"

bash workflow/run_smorfs_pipeline.sh --run-metaeuk ${METAEUK_RUN_ARGS} --cpus "${CPUS}"
