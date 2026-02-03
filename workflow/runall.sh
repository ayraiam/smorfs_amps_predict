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
CONSTRAINT=""

# Optional smORF prediction (OFF by default)
RUN_SMORFS=0
SMORFS_SAMPLE_ID=""        # optional: run only one SampleID

# Optional downstream analysis (MetaFlye flye.log -> table + plots)
RUN_DOWNSTREAM=0
METRICS_ENV="metaflye_metrics_env"

# Output root (per your request)
RESULTS_DIR="results"

# QC-only by default
RUN_FILTERING=0

# Optional QC checks (OFF by default)
RUN_PORECHOP=0

# Optional assembly (OFF by default)
RUN_METAFlyE=0

# smORFs settings
SMORFS_ENV="smorfs_amps_env"
FUNANNOTATE_DB_DIR=""      # optional; exported to funannotate
SMORFS_CREATE_ENV=0
SMORFS_SAMPLES_FILE="metadata/sample_ids.txt"   # default list (one SampleID per line)

# Step control (default: run QC; MetaFlye only if requested)
RUN_QC=1
RUN_ASSEMBLY=0

# Optional filtering toggles (only used if --run-filtering is set)
DO_ADAPTER_TRIM=1
DO_BARCODE_TRIM=1
DO_DEMUX=0
DO_POLY_TRIM=1
DO_QUAL_LEN_FILTER=1

RUN_DOWNSTREAM_AMPS=0
RUN_DOWNSTREAM_MAP=0

MIN_Q="10"
MIN_LEN="500"
MAX_LEN="0"

BATCH_ID="batch1"

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
  echo "Steps:"
  echo "  --qc-only             Run only QC (default)"
  echo "  --metaflye-only       Run only MetaFlye array (skip QC)"
  echo "  --qc-and-metaflye     Run QC then MetaFlye array"
  echo "  --smorfs-only         Run only smORF discovery on existing MetaFlye assemblies"
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
  echo "Assembly:"
  echo "  --run-metaflye        Run metagenome assembly with Flye --meta (default: OFF)"
  echo
  echo "Batch:"
  echo "  --batch-id STR        Batch identifier for NanoPlot/NanoStat labeling (default: batch1)"
  echo "smORFs:"
  echo "  --smorfs-only         Submit smORFs job on existing assemblies (no QC, no assembly)"
  echo "  --smorfs-sample STR   Run smORFs for ONE SampleID only (e.g., TS-0500)"
  echo "  --smorfs-samples FILE Run smORFs for SampleIDs in FILE (one per line)"
  echo "  --smorfs-create-env   Create smORFs conda env and exit"
  echo "  --smorfs-env STR      Conda env name for smORFs pipeline (default: smorfs_amps_env)"
  echo "  --funannotate-db PATH Funannotate DB dir (export FUNANNOTATE_DB_DIR)"
  echo "Downstream:"
  echo "  --downstream-only          Run downstream analysis only (parse flye.log + boxplots)"
  echo "  --run-downstream           Run downstream analysis in addition to selected steps (CAUTION: run after assembly finishes)"
  echo "  --metrics-env STR          Env name for downstream analysis (default: metaflye_metrics_env)"
  echo "  --constraint STR      Slurm constraint (e.g., skylake_avx512)"

  echo
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
    --constraint) CONSTRAINT="$2"; shift 2 ;;

    --qc-only)
      RUN_QC=1
      RUN_ASSEMBLY=0
      RUN_METAFlyE=0
      shift 1
      ;;
    --metaflye-only)
      RUN_QC=0
      RUN_ASSEMBLY=1
      RUN_METAFlyE=1
      shift 1
      ;;
    --qc-and-metaflye)
      RUN_QC=1
      RUN_ASSEMBLY=1
      RUN_METAFlyE=1
      shift 1
      ;;

    --batch-id) BATCH_ID="$2"; shift 2 ;;

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

    --run-metaflye) RUN_METAFlyE=1; shift 1 ;;

    --run-smorfs) RUN_SMORFS=1; shift 1 ;;

    -h|--help) usage ;;

    # ---- smORFs options ----
    --smorfs-only)
      RUN_QC=0
      RUN_ASSEMBLY=0
      RUN_METAFlyE=0
      RUN_SMORFS=1
      shift 1
      ;;
    --smorfs-sample) SMORFS_SAMPLE_ID="$2"; shift 2 ;;
    --smorfs-samples) SMORFS_SAMPLES_FILE="$2"; shift 2 ;;
    --smorfs-create-env) SMORFS_CREATE_ENV=1; shift 1 ;;
    --smorfs-env) SMORFS_ENV="$2"; shift 2 ;;
    --funannotate-db) FUNANNOTATE_DB_DIR="$2"; shift 2 ;;
    --run-downstream) RUN_DOWNSTREAM=1; shift 1 ;;
    --downstream-only)
      RUN_QC=0
      RUN_ASSEMBLY=0
      RUN_METAFlyE=0
      RUN_SMORFS=0
      RUN_DOWNSTREAM=1
      shift 1
      ;;
    --metrics-env) METRICS_ENV="$2"; shift 2 ;;
    --run-downstream-amps) RUN_DOWNSTREAM_AMPS=1; shift 1 ;;
    --run-downstream-map) RUN_DOWNSTREAM_MAP=1; shift 1 ;;
    --downstream-full) RUN_DOWNSTREAM=1; RUN_DOWNSTREAM_AMPS=1; RUN_DOWNSTREAM_MAP=1; shift 1 ;;

    *) echo "Unknown argument: $1"; usage ;;
  esac
done

# If user asked for MetaFlye but didn't specify steps, run QC + assembly
if [[ "${RUN_METAFlyE}" -eq 1 && "${RUN_ASSEMBLY}" -eq 0 ]]; then
  RUN_ASSEMBLY=1
fi

mkdir -p logs metadata

TS=$(date +%Y%m%d_%H%M%S)
OUT_LOG="logs/qc_${TS}.out"
ERR_LOG="logs/qc_${TS}.err"
CMD_LOG="logs/command_${TS}.txt"
MF_OUT_LOG="logs/metaflye_submit_${TS}.out"
MF_ERR_LOG="logs/metaflye_submit_${TS}.err"
SMORFS_OUT_LOG="logs/smorfs_submit_${TS}.out"
SMORFS_ERR_LOG="logs/smorfs_submit_${TS}.err"
DS_OUT_LOG="logs/downstream_${TS}.out"
DS_ERR_LOG="logs/downstream_${TS}.err"

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
  echo "  BATCH_ID       : ${BATCH_ID}"
  echo "  RESULTS_DIR     : ${RESULTS_DIR}"
  echo "  RUN_FILTERING   : ${RUN_FILTERING}"
  echo "  RUN_PORECHOP    : ${RUN_PORECHOP}"
  echo "  RUN_QC          : ${RUN_QC}"
  echo "  RUN_ASSEMBLY    : ${RUN_ASSEMBLY}"
  echo "  RUN_METAFlyE    : ${RUN_METAFlyE}"
  echo "  PARTITION/TIME  : ${PARTITION} / ${TIME}"
  echo "  CPUS/MEM        : ${CPUS} / ${MEM}"
  echo "  RUN_SMORFS      : ${RUN_SMORFS}"
  echo "  SMORFS_ENV      : ${SMORFS_ENV}"
  echo "  FUNANNOTATE_DB  : ${FUNANNOTATE_DB_DIR:-<unset>}"
  echo "  RUN_DOWNSTREAM : ${RUN_DOWNSTREAM}"
  echo "  METRICS_ENV    : ${METRICS_ENV}"
  echo "============================================"
  echo
} | tee -a "$OUT_LOG" "$ERR_LOG" "$CMD_LOG"

export OMP_NUM_THREADS="$CPUS"
export MKL_NUM_THREADS="$CPUS"
export NUMEXPR_NUM_THREADS="$CPUS"
export PIPELINE_INVOCATION="$CMDLINE"

if [[ "${SMORFS_CREATE_ENV}" -eq 1 ]]; then
  echo ">>> Creating smORFs env via workflow/run_smorfs_pipeline.sh --create-env" | tee -a "$OUT_LOG" "$CMD_LOG"
  /bin/bash workflow/run_smorfs_pipeline.sh --create-env \
    >>"$OUT_LOG" 2>>"$ERR_LOG"
  echo ">>> Env creation complete. Exiting as requested (--smorfs-create-env)." | tee -a "$OUT_LOG" "$CMD_LOG"
  exit 0
fi

if [[ "${RUN_QC}" -eq 1 ]]; then
  srun \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME" \
    --chdir="$WDIR" \
    --export=ALL,THREADS="$CPUS",RESULTS_DIR="$RESULTS_DIR",BATCH_ID="$BATCH_ID",RUN_FILTERING="$RUN_FILTERING",RUN_PORECHOP="$RUN_PORECHOP",RUN_METAFlyE="$RUN_METAFlyE",DO_ADAPTER_TRIM="$DO_ADAPTER_TRIM",DO_BARCODE_TRIM="$DO_BARCODE_TRIM",DO_DEMUX="$DO_DEMUX",DO_POLY_TRIM="$DO_POLY_TRIM",DO_QUAL_LEN_FILTER="$DO_QUAL_LEN_FILTER",MIN_Q="$MIN_Q",MIN_LEN="$MIN_LEN",MAX_LEN="$MAX_LEN",PIPELINE_INVOCATION="$PIPELINE_INVOCATION" \
    /bin/bash workflow/run_libsQC.sh \
    >>"$OUT_LOG" \
    2>>"$ERR_LOG"
else
  echo ">>> Skipping QC step (RUN_QC=0)" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_ASSEMBLY}" -eq 1 && "${RUN_METAFlyE}" -eq 1 ]]; then
  echo ">>> Submitting MetaFlye Slurm array (co-assembly per SampleID) ..."
  echo ">>> MetaFlye submit logs:"
  echo "  ${MF_OUT_LOG}"
  echo "  ${MF_ERR_LOG}"

  {
    echo "============================================"
    echo "MetaFlye submission (from runall.sh)"
    echo "--------------------------------------------"
    echo "Timestamp : $(date --iso-8601=seconds)"
    echo "Host      : $(hostname)"
    echo "PWD       : $(pwd)"
    echo
    echo "Invocation:"
    echo "  ${CMDLINE}"
    echo
    echo "Env passed to submit script:"
    echo "  PARTITION    : ${PARTITION}"
    echo "  TIME         : ${TIME}"
    echo "  CPUS         : ${CPUS}"
    echo "  MEM          : ${MEM}"
    echo "  WDIR         : ${WDIR}"
    echo "  RESULTS_DIR  : ${RESULTS_DIR}"
    echo "  FASTQ_DIR    : data"
    echo "  METADATA_MAP : metadata/metagenome_files.txt"
    echo "============================================"
    echo
  } >>"${MF_OUT_LOG}" 2>>"${MF_ERR_LOG}"

  PARTITION="$PARTITION" TIME="$TIME" CPUS="$CPUS" MEM="$MEM" WDIR="$WDIR" \
  RESULTS_DIR="$RESULTS_DIR" FASTQ_DIR="data" METADATA_MAP="metadata/metagenome_files.txt" \
  SAMPLE_COL="SampleID" FASTQ_COL="FASTQ_Filename" \
  bash workflow/submit_metaflye_array.sh \
    >>"${MF_OUT_LOG}" \
    2>>"${MF_ERR_LOG}"

fi

SMORFS_JOB_ID=""

if [[ "${RUN_SMORFS}" -eq 1 ]]; then
  echo ">>> Submitting smORFs job (assumes MetaFlye assemblies already exist) ..."
  echo ">>> smORFs submit logs:"
  echo "  ${SMORFS_OUT_LOG}"
  echo "  ${SMORFS_ERR_LOG}"

  # Decide whether to run a single sample or a samples file
  SMORFS_RUN_ARGS="--samples-file ${SMORFS_SAMPLES_FILE}"
  if [[ -n "${SMORFS_SAMPLE_ID}" ]]; then
    SMORFS_RUN_ARGS="--sample ${SMORFS_SAMPLE_ID}"
    echo ">>> smORFs will run for ONE SampleID only: ${SMORFS_SAMPLE_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  else
    echo ">>> smORFs will run for SampleIDs in: ${SMORFS_SAMPLES_FILE}" | tee -a "$OUT_LOG" "$CMD_LOG"
  fi

  SMORFS_JOB_ID=$(sbatch \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --export=ALL,CPUS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",FUNANNOTATE_DB_DIR="${FUNANNOTATE_DB_DIR}",SMORFS_ENV="${SMORFS_ENV}",SMORFS_RUN_ARGS="${SMORFS_RUN_ARGS}" \
    --output="${SMORFS_OUT_LOG}" \
    --error="${SMORFS_ERR_LOG}" \
    workflow/smorfs_job.sh
  )

  echo ">>> smORFs submitted as job ${SMORFS_JOB_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_DOWNSTREAM}" -eq 1 ]]; then
  echo ">>> Running downstream analysis (flye.log -> metrics + boxplots) ..." | tee -a "$OUT_LOG" "$CMD_LOG"
  echo ">>> Downstream logs:"
  echo "  ${DS_OUT_LOG}"
  echo "  ${DS_ERR_LOG}"

  # Build downstream CLI arguments
  DS_ARGS=( --results-dir "$RESULTS_DIR" --metrics-env "$METRICS_ENV" )

  if [[ "${RUN_DOWNSTREAM_AMPS}" -eq 1 ]]; then
    DS_ARGS+=( --run-amps )
  fi

  if [[ "${RUN_DOWNSTREAM_MAP}" -eq 1 ]]; then
    DS_ARGS+=( --run-map )
  fi

  SRUN_CONSTRAINT=()
  if [[ -n "${CONSTRAINT}" ]]; then
    SRUN_CONSTRAINT=( --constraint="${CONSTRAINT}" )
  fi

  srun \
    "${SRUN_CONSTRAINT[@]}" \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME" \
    --chdir="$WDIR" \
    --export=ALL,RESULTS_DIR="$RESULTS_DIR",METRICS_ENV="$METRICS_ENV",CPUS="$CPUS" \
    /bin/bash workflow/downstream_analysis.sh \
      "${DS_ARGS[@]}" \
    >>"$DS_OUT_LOG" 2>>"$DS_ERR_LOG"
fi

echo ">>> Pipeline finished."
echo "Logs:"
echo "  $OUT_LOG"
echo "  $ERR_LOG"
echo "Command record:"
echo "  $CMD_LOG"
echo "MetaFlye submit logs:"
echo "  $MF_OUT_LOG"
echo "  $MF_ERR_LOG"

