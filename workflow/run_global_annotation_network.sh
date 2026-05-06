#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR="${RESULTS_DIR:-results}"
PARTITION="${PARTITION:-short}"
TIME="${TIME:-04:00:00}"
CPUS="${CPUS:-8}"
MEM="${MEM:-32G}"
WDIR="${WDIR:-$PWD}"

ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ENV_NAME="${ENV_NAME:-global_annot_network_env}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"

ANNOTATION_FILE="${RESULTS_DIR}/annotation/global_cds/global_rep_cds_annotation.tsv"
COUNT_MATRIX="${RESULTS_DIR}/differential_abundance/aldex2/aldex2_counts_matrix.tsv"
METADATA_FILE="${RESULTS_DIR}/differential_abundance/aldex2/aldex2_metadata_used.tsv"
OUTDIR="${RESULTS_DIR}/annotation_network/global_cross_layer_network"

MIN_EDGE_SUPPORT="${MIN_EDGE_SUPPORT:-2}"
MIN_EDGE_SUM_CPM="${MIN_EDGE_SUM_CPM:-0}"
TOP_EDGE_PERCENTILE="${TOP_EDGE_PERCENTILE:-0}"
NETWORK_STEP="${NETWORK_STEP:-all}"

CREATE_ENV=0
RUN_ANALYSIS=1

die(){ echo "ERROR: $*" >&2; exit 1; }
msg(){ echo "[$(date +'%F %T')] $*" >&2; }

usage(){
  cat <<EOF
Usage:
  bash workflow/run_global_annotation_network.sh [options]

Options:
  --create-env                  Create conda env and exit
  --env-name NAME               Conda env name (default: global_annot_network_env)
  --results-dir DIR             Results directory (default: results)
  --annotation-file PATH         Annotation TSV
  --count-matrix PATH            ALDEx2 count matrix TSV
  --metadata-file PATH           Metadata TSV
  --outdir PATH                  Output directory
  --partition STR                Slurm partition
  --time HH:MM:SS                Slurm time
  --cpus INT                     CPUs
  --mem STR                      Memory
  --wd PATH                      Working directory
  --min-edge-support INT         Minimum CDS support per edge (default: 2)
  --min-edge-sum-cpm FLOAT       Minimum summed global CPM per edge (default: 0)
  --top-edge-percentile FLOAT    Optional edge abundance percentile cutoff, e.g. 75
  --step STR                    Step to run: all, 1-10, or named step
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --create-env) CREATE_ENV=1; RUN_ANALYSIS=0; shift ;;
    --env-name) ENV_NAME="$2"; ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --annotation-file) ANNOTATION_FILE="$2"; shift 2 ;;
    --count-matrix) COUNT_MATRIX="$2"; shift 2 ;;
    --metadata-file) METADATA_FILE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --wd) WDIR="$2"; shift 2 ;;
    --min-edge-support) MIN_EDGE_SUPPORT="$2"; shift 2 ;;
    --min-edge-sum-cpm) MIN_EDGE_SUM_CPM="$2"; shift 2 ;;
    --top-edge-percentile) TOP_EDGE_PERCENTILE="$2"; shift 2 ;;
    --step) NETWORK_STEP="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

init_conda() {
  if command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi

  for guess in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "${guess}/etc/profile.d/conda.sh" ]]; then
      source "${guess}/etc/profile.d/conda.sh"
      return 0
    fi
  done

  die "conda not found"
}

create_env() {
  init_conda
  mkdir -p "${ENV_PREFIX_DIR}"

  if [[ -d "${ENV_PREFIX}" && -f "${ENV_PREFIX}/.env_ready" ]]; then
    msg "Env already exists: ${ENV_PREFIX}"
    return 0
  fi

  msg "Creating env: ${ENV_PREFIX}"

  local installer="conda"
  if command -v mamba >/dev/null 2>&1; then
    installer="mamba"
  fi

  rm -rf "${ENV_PREFIX}" 2>/dev/null || true

  "${installer}" create -y -p "${ENV_PREFIX}" -c conda-forge -c bioconda \
    python=3.11 pandas numpy networkx python-louvain

  touch "${ENV_PREFIX}/.env_ready"

  msg "Env created: ${ENV_PREFIX}"
}

run_network() {
  [[ -s "${ANNOTATION_FILE}" ]] || die "Missing annotation file: ${ANNOTATION_FILE}"
  [[ -s "${COUNT_MATRIX}" ]] || die "Missing count matrix: ${COUNT_MATRIX}"

  mkdir -p "${OUTDIR}" logs

  msg "Submitting global annotation network job"

  sbatch \
    --job-name=global_annot_network \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --output="logs/global_annot_network_%j.out" \
    --error="logs/global_annot_network_%j.err" \
    --export=ALL,ENV_PREFIX="${ENV_PREFIX}",ANNOTATION_FILE="${ANNOTATION_FILE}",COUNT_MATRIX="${COUNT_MATRIX}",METADATA_FILE="${METADATA_FILE}",OUTDIR="${OUTDIR}",NETWORK_STEP="${NETWORK_STEP}",MIN_EDGE_SUPPORT="${MIN_EDGE_SUPPORT}",MIN_EDGE_SUM_CPM="${MIN_EDGE_SUM_CPM}",TOP_EDGE_PERCENTILE="${TOP_EDGE_PERCENTILE}" \
    --wrap='
set -euo pipefail
export PATH="${ENV_PREFIX}/bin:${PATH}"
hash -r
python workflow/build_global_cross_layer_network.py \
  --step "${NETWORK_STEP}" \
  --annotation-file "${ANNOTATION_FILE}" \
  --count-matrix "${COUNT_MATRIX}" \
  --metadata-file "${METADATA_FILE}" \
  --outdir "${OUTDIR}" \
  --min-edge-support "${MIN_EDGE_SUPPORT}" \
  --min-edge-sum-cpm "${MIN_EDGE_SUM_CPM}" \
  --top-edge-percentile "${TOP_EDGE_PERCENTILE}"
'
}

if [[ "${CREATE_ENV}" -eq 1 ]]; then
  create_env
fi

if [[ "${RUN_ANALYSIS}" -eq 1 ]]; then
  create_env
  run_network
fi
