#!/usr/bin/env bash
set -euo pipefail

# conda init
source "$(conda info --base)/etc/profile.d/conda.sh"

# Ensure env + tools
bash workflow/run_smorfs_pipeline.sh --ensure-env

# Run
# shellcheck disable=SC2086
bash workflow/run_smorfs_pipeline.sh --run ${SMORFS_RUN_ARGS} --cpus "${CPUS}"
