#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_libsQC.sh
# Purpose: QC + (optional) cleanup for ONT shotgun metagenomes (FASTQ in data/)
# Default: QC-only (NO trimming/filtering)
# Outputs (per your request):
#   results/qc_pre_filt/
# ==========================================================
set -euo pipefail

if ! command -v conda >/dev/null 2>&1; then
  for CAND in "$HOME/mambaforge" "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [ -f "$CAND/etc/profile.d/conda.sh" ]; then
      source "$CAND/etc/profile.d/conda.sh"
      break
    fi
  done
fi

THREADS="${THREADS:-8}"
RESULTS_DIR="${RESULTS_DIR:-results}"
RESULTS="${RESULTS_DIR}"
BATCH_ID="${BATCH_ID:-batch1}"

# Requested pre-QC directory
QC_PRE_DIR="${RESULTS}/qc_pre_filt"
QC_POST_DIR="${RESULTS}/qc_post_filt"

RUN_FILTERING="${RUN_FILTERING:-0}"

# Optional QC checks
RUN_PORECHOP="${RUN_PORECHOP:-0}"

DO_ADAPTER_TRIM="${DO_ADAPTER_TRIM:-1}"
DO_BARCODE_TRIM="${DO_BARCODE_TRIM:-1}"
DO_DEMUX="${DO_DEMUX:-0}"
DO_POLY_TRIM="${DO_POLY_TRIM:-1}"
DO_QUAL_LEN_FILTER="${DO_QUAL_LEN_FILTER:-1}"

MIN_Q="${MIN_Q:-10}"
MIN_LEN="${MIN_LEN:-500}"
MAX_LEN="${MAX_LEN:-0}"

echo "============================================"
echo "QC job started (run_libsQC.sh)"
echo "--------------------------------------------"
echo "Timestamp   : $(date --iso-8601=seconds)"
echo "Host        : $(hostname)"
echo "PWD         : $(pwd)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
echo
echo "Pipeline invocation:"
echo "  ${PIPELINE_INVOCATION:-<not provided>}"
echo
echo "Output roots:"
echo "  RESULTS_DIR : ${RESULTS_DIR}"
echo "  QC_PRE_DIR  : ${QC_PRE_DIR}"
echo "Mode:"
echo "  RUN_FILTERING: ${RUN_FILTERING}"
echo "QC checks:"
echo "  RUN_PORECHOP : ${RUN_PORECHOP}"
echo "Batch:"
echo "  BATCH_ID      : ${BATCH_ID}"
echo "============================================"
echo

ensure_channels() {
  conda config --remove-key channels >/dev/null 2>&1 || true
  conda config --add channels conda-forge
  conda config --add channels bioconda
  conda config --add channels defaults
  conda config --set channel_priority strict
  conda clean --index-cache -y >/dev/null 2>&1 || true
}

create_env_metaQC() {
  local ENV_NAME="metaQC"
  if conda env list | grep -qE "^${ENV_NAME}\s"; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${ENV_NAME}"
    return 0
  fi

  ensure_channels
  local SOLVER="mamba"
  command -v mamba >/dev/null 2>&1 || SOLVER="conda"

  set +e
  ${SOLVER} create -n "${ENV_NAME}" \
    -c conda-forge -c bioconda \
    python=3.11 \
    fastqc=0.12.1 multiqc=1.21 \
    "seqkit>=2.6" \
    nanostat nanoplot nanofilt \
    "cutadapt>=4.5" \
    porechop pigz \
    -y
  status=$?
  set -e

  if [ $status -ne 0 ]; then
    conda config --set channel_priority flexible
    ${SOLVER} create -n "${ENV_NAME}" \
      -c conda-forge -c bioconda \
      python=3.11 \
      fastqc=0.12.1 multiqc=1.21 \
      "seqkit>=2.6" \
      nanostat nanoplot nanofilt \
      "cutadapt>=4.5" \
      porechop pigz \
      -y
    conda config --set channel_priority strict
  fi

  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${ENV_NAME}"
}

export_env() {
  mkdir -p envs
  conda env export --name metaQC > envs/metaQC.yml
}

gather_fastq_files() {
  shopt -s nullglob
  FASTQ_FILES=( data/*.fastq.gz data/*.fq.gz data/*.fastq data/*.fq )
  shopt -u nullglob
  if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "!!! No FASTQ files found in data/"
    exit 1
  fi
  echo ">>> Found ${#FASTQ_FILES[@]} FASTQ files in data/."
}

run_fastqc_all() {
  local OUTDIR="$1"
  mkdir -p "${OUTDIR}"

  # Limit Java heap used by FastQC (tune as needed)
  export _JAVA_OPTIONS="-Xmx2g"

  fastqc -t "${THREADS}" -o "${OUTDIR}" "${FASTQ_FILES[@]}"
}

run_multiqc() {
  local INDIR="$1"
  local OUTDIR="$2"
  mkdir -p "${OUTDIR}"
  multiqc -o "${OUTDIR}" "${INDIR}"
}

quick_len_qual_overview() {
  local OUTROOT="$1"
  local LABEL="${2:-raw}"

  # Batch-aware label/path
  local TAG="${BATCH_ID}_${LABEL}"

  mkdir -p "${OUTROOT}/nanoplot/${TAG}" "${OUTROOT}/nanostat/${TAG}"

  for f in "${FASTQ_FILES[@]}"; do
    local base="$(basename "$f")"
    base="${base%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"

    NanoStat --fastq "$f" --threads "${THREADS}" \
      --outdir "${OUTROOT}/nanostat/${TAG}/${base}.stat" \
      --name "${base}" >/dev/null 2>&1 || true
  done

  NanoPlot --threads "${THREADS}" --fastq "${FASTQ_FILES[@]}" \
    --N50 --loglength --plots hex dot kde --tsv_stats --raw \
    -o "${OUTROOT}/nanoplot/${TAG}" \
    1>"${OUTROOT}/nanoplot/${TAG}/NanoPlot.stdout.log" \
    2>"${OUTROOT}/nanoplot/${TAG}/NanoPlot.stderr.log" || true
}

make_fastq_summary() {
  local OUTROOT="$1"
  local LABEL="${2:-raw}"
  mkdir -p "${OUTROOT}/summary"

  local TAG="${BATCH_ID}_${LABEL}"
  seqkit stats -a -T "${FASTQ_FILES[@]}" > "${OUTROOT}/summary/seqkit_stats_${TAG}.tsv"
}

porechop_check() {
  local OUTROOT="$1"
  mkdir -p "${OUTROOT}/adapter_barcode_checks"
  for f in "${FASTQ_FILES[@]}"; do
    local base="$(basename "$f")"
    base="${base%.fastq.gz}"; base="${base%.fq.gz}"; base="${base%.fastq}"; base="${base%.fq}"
    porechop --check_reads 5000 -i "$f" \
      > "${OUTROOT}/adapter_barcode_checks/${base}.porechop_check.txt" 2>&1 || true
  done
}

time_function() {
  local fn="$1"
  local start=$(date +%s)
  echo ">>> Running $fn ..."
  set +e; $fn; local status=$?; set -e
  local end=$(date +%s); local dur=$((end - start))
  mkdir -p logs
  echo -e "${fn}\t${dur}" >> logs/.timing.tsv
  echo ">>> $fn completed in ${dur}s"
  return $status
}

# ==========================================================
# Pipeline functions
# ==========================================================
pre_qc() {
  echo ">>> === PRE-QC (read-only) ==="
  mkdir -p "${RESULTS}" "${QC_PRE_DIR}"

  time_function "run_fastqc_all ${QC_PRE_DIR}/fastqc"
  time_function "run_multiqc ${QC_PRE_DIR}/fastqc ${QC_PRE_DIR}/multiqc"
  time_function "quick_len_qual_overview ${QC_PRE_DIR} raw"
  time_function "make_fastq_summary ${QC_PRE_DIR} raw"

  if [[ "${RUN_PORECHOP}" -eq 1 ]]; then
    time_function "porechop_check ${QC_PRE_DIR}"
  else
    echo ">>> Skipping porechop_check (RUN_PORECHOP=0)."
  fi
}

run_filtering() {
  echo ">>> === FILTERING/CLEANUP (optional) ==="
  echo ">>> (Not implemented in this QC-only adjustment; enable later if needed.)"
  return 0
}

post_qc() {
  echo ">>> === POST-QC (optional) ==="
  mkdir -p "${QC_POST_DIR}"
  time_function "run_fastqc_all ${QC_POST_DIR}/fastqc"
  time_function "run_multiqc ${QC_POST_DIR}/fastqc ${QC_POST_DIR}/multiqc"
  time_function "quick_len_qual_overview ${QC_POST_DIR} clean"
  time_function "make_fastq_summary ${QC_POST_DIR} clean"
}

# ==========================================================
# Main flow
# ==========================================================
START_TIME=$(date +%s)
mkdir -p logs
: > logs/.timing.tsv

time_function create_env_metaQC
time_function export_env
time_function gather_fastq_files

pre_qc

if [[ "${RUN_FILTERING}" -eq 1 ]]; then
  run_filtering
  post_qc
else
  echo ">>> RUN_FILTERING=0: QC-only run completed. No files were modified."
fi

echo ">>> Done."
echo ">>> Pre-QC outputs: ${QC_PRE_DIR}"
