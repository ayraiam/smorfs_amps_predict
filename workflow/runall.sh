#!/usr/bin/env bash
# ==========================================================
# Script: workflow/runall.sh
# Purpose: Submit ONT metagenome QC via Slurm (QC-only by default).
# Input: FASTQ(.gz) in data/
# Output:
#   results/qc_pre_filt/   (always)
#   results/qc_post_filt/  (only if filtering enabled)
#   results/trimmed|polytrim|filtered (only if filtering enabled)
# ==========================================================
set -euo pipefail

ORIG_ARGS=("$@")

PARTITION="short"
TIME="04:00:00"
CPUS="8"
MEM="32G"
WDIR="$PWD"

# Output root (per your request)
RESULTS_DIR="results"

# QC-only by default
RUN_FILTERING=0

# Optional QC checks (OFF by default)
RUN_PORECHOP=0

# Optional filtering toggles (only used if --run-filtering is set)
DO_ADAPTER_TRIM=1
DO_BARCODE_TRIM=1
DO_DEMUX=0
DO_POLY_TRIM=1
DO_QUAL_LEN_FILTER=1

MIN_Q="10"
MIN_LEN="500"
MAX_LEN="0"

usage() {
  echo "Usage: bash workflow/runall.sh [options]"
  echo
  echo "Resources:"
  echo "  --partition STR       (default: short)"
  echo "  --time HH:MM:SS       (default: 04:00:00)"
  echo "  --cpus INT            (default: 8)"
  echo "  --mem STR             (default: 32G)"
  echo "  --wd PATH             (default: current dir)"
  echo
  echo "Mode:"
  echo "  --run-filtering       Enable trimming/filtering step (default: OFF; QC-only)"
  echo
  echo "QC checks:"
  echo "  --run-porechop        Run porechop adapter/barcode check step (default: OFF)"
  echo
  echo "Optional filtering toggles (only used if --run-filtering):"
  echo "  --no-adapter-trim     Skip adapter trimming"
  echo "  --no-barcode-trim     Skip barcode trimming"
  echo "  --demux               Enable demultiplexing"
  echo "  --no-poly-trim        Skip poly-A/T trimming"
  echo "  --no-filter           Skip NanoFilt Q/len filtering"
  echo "  --min-q INT           Mean read Q cutoff (default: 10)"
  echo "  --min-len INT         Min length bp (default: 500)"
  echo "  --max-len INT         Max length bp (0 disables; default: 0)"
  echo
  echo "Output:"
  echo "  --results-dir PATH    Output root (default: results)"
  echo
  echo "  -h, --help            Show help"
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --wd) WDIR="$2"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;

    --run-filtering) RUN_FILTERING=1; shift 1 ;;

    --run-porechop) RUN_PORECHOP=1; shift 1 ;;

    --no-adapter-trim) DO_ADAPTER_TRIM=0; shift 1 ;;
    --no-barcode-trim) DO_BARCODE_TRIM=0; shift 1 ;;
    --demux) DO_DEMUX=1; shift 1 ;;
    --no-poly-trim) DO_POLY_TRIM=0; shift 1 ;;
    --no-filter) DO_QUAL_LEN_FILTER=0; shift 1 ;;
    --min-q) MIN_Q="$2"; shift 2 ;;
    --min-len) MIN_LEN="$2"; shift 2 ;;
    --max-len) MAX_LEN="$2"; shift 2 ;;

    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

mkdir -p logs metadata

TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="logs/qc_${TS}.out"
ERR_LOG="logs/qc_${TS}.err"
CMD_LOG="logs/command_${TS}.txt"

CMDLINE="$(printf "%q " "$0" "${ORIG_ARGS[@]}")"
RUN_TS="$(date --iso-8601=seconds)"
RUN_HOST="$(hostname)"
RUN_PWD="$(pwd)"

{
  echo "============================================"
  echo "Pipeline invocation (from runall.sh)"
  echo "--------------------------------------------"
  echo "Timestamp : ${RUN_TS}"
  echo "Host      : ${RUN_HOST}"
  echo "PWD       : ${RUN_PWD}"
  echo
  echo "Command:"
  echo "  ${CMDLINE}"
  echo
  echo "Resolved settings:"
  echo "  RESULTS_DIR     : ${RESULTS_DIR}"
  echo "  RUN_FILTERING   : ${RUN_FILTERING}"
  echo "  RUN_PORECHOP    : ${RUN_PORECHOP}"
  echo "  PARTITION/TIME  : ${PARTITION} / ${TIME}"
  echo "  CPUS/MEM        : ${CPUS} / ${MEM}"
  echo "============================================"
  echo
} | tee -a "$OUT_LOG" "$ERR_LOG" "$CMD_LOG"

export OMP_NUM_THREADS="$CPUS"
export MKL_NUM_THREADS="$CPUS"
export NUMEXPR_NUM_THREADS="$CPUS"
export PIPELINE_INVOCATION="$CMDLINE"

srun \
  --partition="$PARTITION" \
  --nodes=1 \
  --ntasks=1 \
  --cpus-per-task="$CPUS" \
  --mem="$MEM" \
  --time="$TIME" \
  --chdir="$WDIR" \
  --export=ALL,THREADS="$CPUS",RESULTS_DIR="$RESULTS_DIR",RUN_FILTERING="$RUN_FILTERING",RUN_PORECHOP="$RUN_PORECHOP",DO_ADAPTER_TRIM="$DO_ADAPTER_TRIM",DO_BARCODE_TRIM="$DO_BARCODE_TRIM",DO_DEMUX="$DO_DEMUX",DO_POLY_TRIM="$DO_POLY_TRIM",DO_QUAL_LEN_FILTER="$DO_QUAL_LEN_FILTER",MIN_Q="$MIN_Q",MIN_LEN="$MIN_LEN",MAX_LEN="$MAX_LEN",PIPELINE_INVOCATION="$PIPELINE_INVOCATION" \
  /bin/bash workflow/run_libsQC.sh \
  >>"$OUT_LOG" \
  2>>"$ERR_LOG"

echo ">>> Pipeline finished."
echo "Logs:"
echo "  $OUT_LOG"
echo "  $ERR_LOG"
echo "Command record:"
echo "  $CMD_LOG"
