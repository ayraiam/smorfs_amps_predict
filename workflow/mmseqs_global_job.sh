#!/usr/bin/env bash
set -euo pipefail

# --------------------------------------------------
# Ensure refine_annot_smorf_bacs_env exists
# and contains mmseqs2
# --------------------------------------------------

ENV_PREFIX="envs/refine_annot_smorf_bacs_env"

init_conda() {
  if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi
  for guess in "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "${guess}/etc/profile.d/conda.sh" ]]; then
      source "${guess}/etc/profile.d/conda.sh"
      return 0
    fi
  done
  echo "ERROR: conda not found."
  exit 1
}

init_conda

if [[ ! -d "${ENV_PREFIX}" ]]; then
  echo "[INFO] Creating refine_annot_smorf_bacs_env (with mmseqs2)..."
  conda create -y -p "${ENV_PREFIX}" --override-channels \
    -c conda-forge -c bioconda \
    python=3.10 \
    pandas \
    numpy \
    pyarrow \
    biopython \
    gffutils \
    pyranges \
    intervaltree \
    seqkit \
    mmseqs2
else
  echo "[INFO] Using existing refine_annot_smorf_bacs_env"
fi

conda activate "${ENV_PREFIX}"

if ! command -v mmseqs >/dev/null 2>&1; then
  echo "[INFO] mmseqs2 not found inside env. Installing..."
  conda install -y -c conda-forge -c bioconda mmseqs2
fi

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
