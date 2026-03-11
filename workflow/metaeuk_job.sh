#!/usr/bin/env bash
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

# Reuse the same env builder/activation logic from run_smorfs_pipeline.sh
bash workflow/run_smorfs_pipeline.sh --run-metaeuk ${METAEUK_RUN_ARGS} --cpus "${CPUS}"
