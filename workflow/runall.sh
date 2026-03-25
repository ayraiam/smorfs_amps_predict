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

# Output root
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

# ---------------------------
# refine_annot_smorf_bacs (bacterial refinement)
# ---------------------------
RUN_REFINE_BACS=0
REFINE_BACS_CREATE_ENV=0
REFINE_BACS_ENV="refine_annot_smorf_bacs_env"
REFINE_BACS_SAMPLES_FILE="metadata/sample_ids.txt"   # default list (same style)
REFINE_BACS_SAMPLE_ID=""                              # optional single sample
REFINE_BACS_INPUT_TSV=""                              # optional override (advanced)
REFINE_BACS_JOB_ID=""
REFINE_CLUSTER_ONLY=0
REFINE_EUKS_CLUSTER_ONLY=0

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

RUN_MACREL_ATTACH=0
PREDICTED_SMORFS_TSV=""
MACREL_ID_COL="feature_id"
MACREL_SEQ_COL=""
MACREL_ATTACH_OUT=""

MIN_Q="10"
MIN_LEN="500"
MAX_LEN="0"

BATCH_ID="batch1"

# Global MetaFlye assembly
GLOBAL_ASSEMBLY=0
GLOBAL_ID="GLOBAL"

RUN_MMSEQS_GLOBAL=0
MIN_SEQ_ID="0.95"
MMSEQS_COV="0.8"
MMSEQS_COV_MODE="1"

# Optional MetaEuk annotation for fungal/euk contigs (OFF by default)
RUN_METAEUK=0
METAEUK_ENV="${SMORFS_ENV:-smorfs_amps_env}"   # same env as smORFs pipeline
METAEUK_SAMPLES_FILE="metadata/sample_ids.txt"
METAEUK_SAMPLE_ID=""
METAEUK_DB=""
METAEUK_JOB_ID=""

#Refining Euk annotations
RUN_REFINE_EUKS=0
REFINE_EUKS_CREATE_ENV=0
REFINE_EUKS_ENV="refine_annot_smorf_euks_env"
REFINE_EUKS_SAMPLES_FILE="metadata/sample_ids.txt"
REFINE_EUKS_SAMPLE_ID=""
REFINE_EUKS_INPUT_TSV=""
REFINE_EUKS_JOB_ID=""
REFINE_EUKS_STEP1=1
REFINE_EUKS_STEP2=1
REFINE_EUKS_STEP3=1

# Scratch root for sample-level smORFs outputs
SMORFS_WORK_ROOT="${SMORFS_WORK_ROOT:-/scratch/t.sousa/data_used/smorfs}"

# Optional global CDS abundance mapping (OFF by default)
RUN_MAP_GLOBAL_CDS=0
MAP_GLOBAL_CDS_BUILD_REF_ONLY=0
MAP_GLOBAL_CDS_ONLY=0
MAP_GLOBAL_CDS_SAMPLE_ID=""
ABUND_ENV_NAME="smorf_abundance_env"

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
  echo "--refine-bacs-only     Run only bacterial refinement (post Macrel/smORFs)"
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
  echo "  --global                 Run a single global co-assembly using all FASTQs in data/"
  echo "  --global-id STR           SampleID name for global assembly (default: GLOBAL)"
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
  echo "MetaEuk (fungi/euk contigs):"
  echo "  --metaeuk-only            Run only MetaEuk annotation on existing fungi_contigs.fasta"
  echo "  --run-metaeuk             Run MetaEuk stage in addition to selected steps"
  echo "  --metaeuk-sample STR      Run MetaEuk for ONE SampleID only"
  echo "  --metaeuk-samples FILE    Run MetaEuk for SampleIDs in FILE (one per line)"
  echo "  --metaeuk-db PATH         Protein reference DB/FASTA for MetaEuk (required)"
  echo "Refine (fungi/euks):"
  echo "  --refine-euks-only         Submit fungal/euk refinement only"
  echo "  --run-refine-euks          Run fungal/euk refine stage in addition to selected steps"
  echo "  --refine-euks-sample STR   Run fungal/euk refine for ONE SampleID only"
  echo "  --refine-euks-samples FILE Run fungal/euk refine for SampleIDs in FILE"
  echo "  --refine-euks-create-env   Create fungal/euk refine env and exit"
  echo "  --refine-euks-env STR      Conda env name (default: refine_annot_smorf_euks_env)"
  echo "  --refine-euks-input-tsv P  Optional override TSV path"
  echo "  --refine-euks-step1-only    Run only fungal/euk Step 1"
  echo "  --refine-euks-step2-only    Run only fungal/euk Step 2"
  echo "  --refine-euks-step3-only    Run only fungal/euk Step 3"
  echo "  --refine-euks-steps STR     Comma-separated fungal/euk steps (e.g. 1,2,3 or 1,3)"
  echo "Clustering (MMseqs2):"
  echo "  --mmseqs-global-only      Build global peptides FASTA + run MMseqs2 clustering"
  echo "  --run-mmseqs-global       Run MMseqs2 global clustering in addition to selected steps"
  echo "  --min-seq-id FLOAT        MMseqs min identity (default: 0.95)"
  echo "  --mmseqs-cov FLOAT        MMseqs coverage -c (default: 0.8)"
  echo "  --mmseqs-cov-mode INT     MMseqs cov-mode (default: 1)"
  echo "Refine (bacteria):"
  echo "--refine-bacs-only         Submit refine job only (no QC, no assembly, no smORFs)"
  echo "--run-refine-bacs          Run refine stage in addition to selected steps"
  echo "--refine-bacs-sample STR   Run refine for ONE SampleID only"
  echo "--refine-bacs-samples FILE Run refine for SampleIDs in FILE (one per line)"
  echo "--refine-bacs-create-env   Create refine env and exit"
  echo "--refine-bacs-env STR      Conda env name (default: refine_annot_smorf_bacs_env)"
  echo "--refine-bacs-input-tsv P  Optional override TSV path (advanced; usually leave empty)"
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

    --global)
      GLOBAL_ASSEMBLY=1
      shift 1
      ;;
    --global-id)
      GLOBAL_ID="$2"
      shift 2
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

    --macrel-attach-only)
      RUN_QC=0
      RUN_ASSEMBLY=0
      RUN_METAFlyE=0
      RUN_SMORFS=0
      RUN_DOWNSTREAM=1
      RUN_MACREL_ATTACH=1
      shift 1
      ;;
    --predicted-smorfs)
      PREDICTED_SMORFS_TSV="$2"
      shift 2
      ;;
    --id-col)
      MACREL_ID_COL="$2"
      shift 2
      ;;
    --seq-col)
      MACREL_SEQ_COL="$2"
      shift 2
      ;;
    --macrel-attach-out)
      MACREL_ATTACH_OUT="$2"
      shift 2
      ;;
			# ---- refine bacs options ----
      --refine-bacs-only)
        RUN_QC=0
        RUN_ASSEMBLY=0
        RUN_METAFlyE=0
        RUN_SMORFS=0
        RUN_DOWNSTREAM=0
        RUN_REFINE_BACS=1
        shift 1
        ;;

      --run-refine-bacs) RUN_REFINE_BACS=1; shift 1 ;;
      --refine-bacs-sample) REFINE_BACS_SAMPLE_ID="$2"; shift 2 ;;
      --refine-bacs-samples) REFINE_BACS_SAMPLES_FILE="$2"; shift 2 ;;
      --refine-bacs-create-env) REFINE_BACS_CREATE_ENV=1; shift 1 ;;
      --refine-bacs-env) REFINE_BACS_ENV="$2"; shift 2 ;;
      --refine-bacs-input-tsv) REFINE_BACS_INPUT_TSV="$2"; shift 2 ;;
      --cluster-only) REFINE_CLUSTER_ONLY=1; shift 1 ;;

      --mmseqs-global-only)
        RUN_QC=0
        RUN_ASSEMBLY=0
        RUN_METAFlyE=0
        RUN_SMORFS=0
        RUN_DOWNSTREAM=0
        RUN_REFINE_BACS=0
        RUN_MMSEQS_GLOBAL=1
        shift 1
        ;;
      --run-mmseqs-global) RUN_MMSEQS_GLOBAL=1; shift 1 ;;
      --min-seq-id) MIN_SEQ_ID="$2"; shift 2 ;;
      --mmseqs-cov) MMSEQS_COV="$2"; shift 2 ;;
      --mmseqs-cov-mode) MMSEQS_COV_MODE="$2"; shift 2 ;;

      # ---- MetaEuk options ----
      --metaeuk-only)
        RUN_QC=0
        RUN_ASSEMBLY=0
        RUN_METAFlyE=0
        RUN_SMORFS=0
        RUN_DOWNSTREAM=0
        RUN_REFINE_BACS=0
        RUN_METAEUK=1
        shift 1
        ;;
      --run-metaeuk) RUN_METAEUK=1; shift 1 ;;
      --metaeuk-sample) METAEUK_SAMPLE_ID="$2"; shift 2 ;;
      --metaeuk-samples) METAEUK_SAMPLES_FILE="$2"; shift 2 ;;
      --metaeuk-db) METAEUK_DB="$2"; shift 2 ;;

      --refine-euks-only)
        RUN_QC=0
        RUN_ASSEMBLY=0
        RUN_METAFlyE=0
        RUN_SMORFS=0
        RUN_DOWNSTREAM=0
        RUN_REFINE_BACS=0
        RUN_REFINE_EUKS=1
        shift 1
        ;;
      --run-refine-euks) RUN_REFINE_EUKS=1; shift 1 ;;
      --refine-euks-sample) REFINE_EUKS_SAMPLE_ID="$2"; shift 2 ;;
      --refine-euks-samples) REFINE_EUKS_SAMPLES_FILE="$2"; shift 2 ;;
      --refine-euks-create-env) REFINE_EUKS_CREATE_ENV=1; shift 1 ;;
      --refine-euks-env) REFINE_EUKS_ENV="$2"; shift 2 ;;
      --refine-euks-input-tsv) REFINE_EUKS_INPUT_TSV="$2"; shift 2 ;;

      --refine-euks-step1-only)
        REFINE_EUKS_STEP1=1
        REFINE_EUKS_STEP2=0
        REFINE_EUKS_STEP3=0
        shift 1
        ;;
      --refine-euks-step2-only)
        REFINE_EUKS_STEP1=0
        REFINE_EUKS_STEP2=1
        REFINE_EUKS_STEP3=0
        shift 1
        ;;
      --refine-euks-step3-only)
        REFINE_EUKS_STEP1=0
        REFINE_EUKS_STEP2=0
        REFINE_EUKS_STEP3=1
        shift 1
        ;;
      --refine-euks-steps)
        REFINE_EUKS_STEP1=0
        REFINE_EUKS_STEP2=0
        REFINE_EUKS_STEP3=0
        IFS=',' read -ra _steps <<< "$2"
        for s in "${_steps[@]}"; do
          [[ "$s" == "1" ]] && REFINE_EUKS_STEP1=1
          [[ "$s" == "2" ]] && REFINE_EUKS_STEP2=1
          [[ "$s" == "3" ]] && REFINE_EUKS_STEP3=1
        done
        shift 2
        ;;

      --refine-euks-cluster-only) REFINE_EUKS_CLUSTER_ONLY=1; shift 1 ;;

      --map-global-cds) RUN_MAP_GLOBAL_CDS=1; shift ;;
      --map-global-cds-build-ref-only) RUN_MAP_GLOBAL_CDS=1; MAP_GLOBAL_CDS_BUILD_REF_ONLY=1; shift ;;
      --map-global-cds-only) RUN_MAP_GLOBAL_CDS=1; MAP_GLOBAL_CDS_ONLY=1; shift ;;
      --map-global-cds-sample-id) RUN_MAP_GLOBAL_CDS=1; MAP_GLOBAL_CDS_SAMPLE_ID="$2"; shift 2 ;;
      --abund-env-name) ABUND_ENV_NAME="$2"; shift 2 ;;

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
  echo "  RUN_METAEUK     : ${RUN_METAEUK}"
  echo "  METAEUK_DB      : ${METAEUK_DB:-<unset>}"
  echo "  SMORFS_WORK_ROOT: ${SMORFS_WORK_ROOT}"
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

if [[ "${REFINE_BACS_CREATE_ENV}" -eq 1 ]]; then
  echo ">>> Creating refine_bacs env via workflow/run_refine_annot_smorf_bacs.sh --create-env" | tee -a "$OUT_LOG" "$CMD_LOG"
  /bin/bash workflow/run_refine_annot_smorf_bacs.sh --create-env \
    --refine-env "${REFINE_BACS_ENV}" \
    >>"$OUT_LOG" 2>>"$ERR_LOG"
  echo ">>> Refine env creation complete. Exiting as requested (--refine-bacs-create-env)." | tee -a "$OUT_LOG" "$CMD_LOG"
  exit 0
fi

if [[ "${REFINE_EUKS_CREATE_ENV}" -eq 1 ]]; then
  echo ">>> Creating refine_euks env via workflow/run_refine_annot_smorf_euks.sh --create-env" | tee -a "$OUT_LOG" "$CMD_LOG"
  /bin/bash workflow/run_refine_annot_smorf_euks.sh --create-env \
    --refine-env "${REFINE_EUKS_ENV}" \
    >>"$OUT_LOG" 2>>"$ERR_LOG"
  echo ">>> Refine euks env creation complete. Exiting as requested (--refine-euks-create-env)." | tee -a "$OUT_LOG" "$CMD_LOG"
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

  if [[ "${GLOBAL_ASSEMBLY}" -eq 1 ]]; then
    echo ">>> MetaFlye GLOBAL co-assembly (SampleID=${GLOBAL_ID})" | tee -a "$MF_OUT_LOG"

    PARTITION="$PARTITION" TIME="$TIME" CPUS="$CPUS" MEM="$MEM" WDIR="$WDIR" \
    RESULTS_DIR="$RESULTS_DIR" FASTQ_DIR="data" \
    bash workflow/submit_metaflye_array.sh \
      --global --global-id "${GLOBAL_ID}" \
      >>"${MF_OUT_LOG}" \
      2>>"${MF_ERR_LOG}"
  else
    PARTITION="$PARTITION" TIME="$TIME" CPUS="$CPUS" MEM="$MEM" WDIR="$WDIR" \
    RESULTS_DIR="$RESULTS_DIR" FASTQ_DIR="data" METADATA_MAP="metadata/metagenome_files.txt" \
    SAMPLE_COL="SampleID" FASTQ_COL="FASTQ_Filename" \
    bash workflow/submit_metaflye_array.sh \
      >>"${MF_OUT_LOG}" \
      2>>"${MF_ERR_LOG}"
  fi

fi

SMORFS_JOB_ID=""
MM_JOB_ID=""

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
    --export=ALL,CPUS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",SMORFS_WORK_ROOT="${SMORFS_WORK_ROOT}",FUNANNOTATE_DB_DIR="${FUNANNOTATE_DB_DIR}",SMORFS_ENV="${SMORFS_ENV}",SMORFS_RUN_ARGS="${SMORFS_RUN_ARGS}" \
    --output="${SMORFS_OUT_LOG}" \
    --error="${SMORFS_ERR_LOG}" \
    workflow/smorfs_job.sh \
    | awk '{print $NF}'
  )

  echo ">>> smORFs submitted as job ${SMORFS_JOB_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_METAEUK}" -eq 1 ]]; then
  [[ -n "${METAEUK_DB}" ]] || { echo "ERROR: --metaeuk-db is required for MetaEuk runs."; exit 1; }

  METAEUK_OUT_LOG="logs/metaeuk_submit_${TS}.out"
  METAEUK_ERR_LOG="logs/metaeuk_submit_${TS}.err"

  echo ">>> Submitting MetaEuk job ..." | tee -a "$OUT_LOG" "$CMD_LOG"

  METAEUK_RUN_ARGS="--samples-file ${METAEUK_SAMPLES_FILE}"
  if [[ -n "${METAEUK_SAMPLE_ID}" ]]; then
    METAEUK_RUN_ARGS="--sample ${METAEUK_SAMPLE_ID}"
    echo ">>> MetaEuk will run for ONE SampleID only: ${METAEUK_SAMPLE_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  else
    echo ">>> MetaEuk will run for SampleIDs in: ${METAEUK_SAMPLES_FILE}" | tee -a "$OUT_LOG" "$CMD_LOG"
  fi

  METAEUK_DEP=()
  if [[ -n "${SMORFS_JOB_ID:-}" ]]; then
    METAEUK_DEP+=( --dependency="afterok:${SMORFS_JOB_ID}" )
  fi

  METAEUK_JOB_ID=$(sbatch \
    "${METAEUK_DEP[@]}" \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --export=ALL,CPUS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",SMORFS_WORK_ROOT="${SMORFS_WORK_ROOT}",SMORFS_ENV="${SMORFS_ENV}",METAEUK_DB="${METAEUK_DB}",METAEUK_RUN_ARGS="${METAEUK_RUN_ARGS}" \
    --output="${METAEUK_OUT_LOG}" \
    --error="${METAEUK_ERR_LOG}" \
    workflow/metaeuk_job.sh \
    | awk '{print $NF}'
  )

  echo ">>> MetaEuk submitted as job ${METAEUK_JOB_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  echo ">>> MetaEuk logs: ${METAEUK_OUT_LOG} / ${METAEUK_ERR_LOG}" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_MMSEQS_GLOBAL}" -eq 1 ]]; then
  MM_OUT_LOG="logs/mmseqs_global_${TS}.out"
  MM_ERR_LOG="logs/mmseqs_global_${TS}.err"

  MM_DEP=()
  if [[ -n "${SMORFS_JOB_ID:-}" ]]; then
    MM_DEP=( --dependency="afterok:${SMORFS_JOB_ID}" )
  fi

  MM_JOB_ID=$(sbatch \
    "${MM_DEP[@]}" \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --export=ALL,RESULTS_DIR="${RESULTS_DIR}",THREADS="${CPUS}",MIN_SEQ_ID="${MIN_SEQ_ID}",COV="${MMSEQS_COV}",COV_MODE="${MMSEQS_COV_MODE}" \
    --output="${MM_OUT_LOG}" \
    --error="${MM_ERR_LOG}" \
    workflow/mmseqs_global_job.sh \
    | awk '{print $NF}'
  )

  echo ">>> MMseqs GLOBAL clustering submitted as job ${MM_JOB_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  echo ">>> MMseqs logs: ${MM_OUT_LOG} / ${MM_ERR_LOG}" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_MAP_GLOBAL_CDS}" -eq 1 ]]; then
  echo "[INFO] Running GLOBAL CDS abundance mapping step"

  map_args=(
    --results-dir "${RESULTS_DIR}"
    --metadata-map "metadata/metagenome_files.txt"
    --partition "${PARTITION}"
    --time "${TIME}"
    --cpus "${CPUS}"
    --mem "${MEM}"
    --env-name "${ABUND_ENV_NAME}"
  )

  if [[ "${MAP_GLOBAL_CDS_BUILD_REF_ONLY}" -eq 1 ]]; then
    map_args+=(--build-ref-only)
  fi

  if [[ "${MAP_GLOBAL_CDS_ONLY}" -eq 1 ]]; then
    map_args+=(--map-only)
  fi

  if [[ -n "${MAP_GLOBAL_CDS_SAMPLE_ID}" ]]; then
    map_args+=(--sample-id "${MAP_GLOBAL_CDS_SAMPLE_ID}")
  fi

  bash workflow/run_map_global_cds_abundance.sh "${map_args[@]}"
fi

if [[ "${RUN_REFINE_BACS}" -eq 1 ]]; then
  echo ">>> Submitting refine_annot_smorf_bacs job ..." | tee -a "$OUT_LOG" "$CMD_LOG"

  REFINE_RUN_ARGS="--samples-file ${REFINE_BACS_SAMPLES_FILE}"

  if [[ -n "${REFINE_BACS_SAMPLE_ID}" ]]; then
    REFINE_RUN_ARGS="--sample ${REFINE_BACS_SAMPLE_ID}"
    echo ">>> refine will run for ONE SampleID only: ${REFINE_BACS_SAMPLE_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  else
    echo ">>> refine will run for SampleIDs in: ${REFINE_BACS_SAMPLES_FILE}" | tee -a "$OUT_LOG" "$CMD_LOG"
  fi

  # Append cluster-only flag if requested (applies to both modes)
  if [[ "${REFINE_CLUSTER_ONLY}" -eq 1 ]]; then
    REFINE_RUN_ARGS="${REFINE_RUN_ARGS} --cluster-only"
    echo ">>> refine running in CLUSTER-ONLY mode (Step 3 only)" | tee -a "$OUT_LOG" "$CMD_LOG"
  fi

  # Optional override TSV (advanced)
  if [[ -n "${REFINE_BACS_INPUT_TSV}" ]]; then
    REFINE_RUN_ARGS="${REFINE_RUN_ARGS} --input-tsv ${REFINE_BACS_INPUT_TSV}"
  fi

  REFINE_OUT_LOG="logs/refine_bacs_submit_${TS}.out"
  REFINE_ERR_LOG="logs/refine_bacs_submit_${TS}.err"

  # If smORFs ran in the same runall.sh invocation, make refine depend on it.
  # SMORFS_JOB_ID will be just the numeric id if you applied the awk fix.
  REFINE_DEP=()
  if [[ -n "${SMORFS_JOB_ID:-}" ]]; then
    REFINE_DEP+=( --dependency="afterok:${SMORFS_JOB_ID}" )
  fi
  if [[ -n "${MM_JOB_ID:-}" ]]; then
    REFINE_DEP+=( --dependency="afterok:${MM_JOB_ID}" )
  fi

  REFINE_JOB_ID=$(sbatch \
    "${REFINE_DEP[@]}" \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --export=ALL,CPUS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",REFINE_BACS_ENV="${REFINE_BACS_ENV}",REFINE_RUN_ARGS="${REFINE_RUN_ARGS}" \
    --output="${REFINE_OUT_LOG}" \
    --error="${REFINE_ERR_LOG}" \
    workflow/refine_bacs_job.sh \
    | awk '{print $NF}'
  )

  echo ">>> refine_annot_smorf_bacs submitted as job ${REFINE_JOB_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  echo ">>> refine logs: ${REFINE_OUT_LOG} / ${REFINE_ERR_LOG}" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_REFINE_EUKS}" -eq 1 ]]; then
  echo ">>> Submitting refine_annot_smorf_euks job ..." | tee -a "$OUT_LOG" "$CMD_LOG"

  REFINE_EUKS_RUN_ARGS="--samples-file ${REFINE_EUKS_SAMPLES_FILE}"

  if [[ -n "${REFINE_EUKS_SAMPLE_ID}" ]]; then
    REFINE_EUKS_RUN_ARGS="--sample ${REFINE_EUKS_SAMPLE_ID}"
    echo ">>> refine_euks will run for ONE SampleID only: ${REFINE_EUKS_SAMPLE_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  else
    echo ">>> refine_euks will run for SampleIDs in: ${REFINE_EUKS_SAMPLES_FILE}" | tee -a "$OUT_LOG" "$CMD_LOG"
  fi

  REFINE_EUKS_RUN_ARGS="${REFINE_EUKS_RUN_ARGS} --run-step1 ${REFINE_EUKS_STEP1} --run-step2 ${REFINE_EUKS_STEP2} --run-step3 ${REFINE_EUKS_STEP3}"

  echo ">>> refine_euks steps: step1=${REFINE_EUKS_STEP1} step2=${REFINE_EUKS_STEP2} step3=${REFINE_EUKS_STEP3}" | tee -a "$OUT_LOG" "$CMD_LOG"

  if [[ "${REFINE_EUKS_CLUSTER_ONLY}" -eq 1 ]]; then
    REFINE_EUKS_RUN_ARGS="${REFINE_EUKS_RUN_ARGS} --cluster-only"
    echo ">>> refine_euks running in CLUSTER-ONLY mode (Step 3 only)" | tee -a "$OUT_LOG" "$CMD_LOG"
  fi

  if [[ -n "${REFINE_EUKS_INPUT_TSV}" ]]; then
    REFINE_EUKS_RUN_ARGS="${REFINE_EUKS_RUN_ARGS} --input-tsv ${REFINE_EUKS_INPUT_TSV}"
  fi

  REFINE_EUKS_OUT_LOG="logs/refine_euks_submit_${TS}.out"
  REFINE_EUKS_ERR_LOG="logs/refine_euks_submit_${TS}.err"

  REFINE_EUKS_DEP=()
  if [[ -n "${METAEUK_JOB_ID:-}" ]]; then
    REFINE_EUKS_DEP+=( --dependency="afterok:${METAEUK_JOB_ID}" )
  fi
  if [[ -n "${MM_JOB_ID:-}" ]]; then
    REFINE_EUKS_DEP+=( --dependency="afterok:${MM_JOB_ID}" )
  fi

  REFINE_EUKS_JOB_ID=$(sbatch \
    "${REFINE_EUKS_DEP[@]}" \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --export=ALL,CPUS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",SMORFS_WORK_ROOT="${SMORFS_WORK_ROOT}",REFINE_EUKS_ENV="${REFINE_EUKS_ENV}",REFINE_EUKS_RUN_ARGS="${REFINE_EUKS_RUN_ARGS}" \
    --output="${REFINE_EUKS_OUT_LOG}" \
    --error="${REFINE_EUKS_ERR_LOG}" \
    workflow/refine_euks_job.sh \
    | awk '{print $NF}'
  )

  echo ">>> refine_annot_smorf_euks submitted as job ${REFINE_EUKS_JOB_ID}" | tee -a "$OUT_LOG" "$CMD_LOG"
  echo ">>> refine_euks logs: ${REFINE_EUKS_OUT_LOG} / ${REFINE_EUKS_ERR_LOG}" | tee -a "$OUT_LOG" "$CMD_LOG"
fi

if [[ "${RUN_DOWNSTREAM}" -eq 1 ]]; then
  echo ">>> Running downstream analysis (flye.log -> metrics + boxplots) ..." | tee -a "$OUT_LOG" "$CMD_LOG"
  echo ">>> Downstream logs:"
  echo "  ${DS_OUT_LOG}"
  echo "  ${DS_ERR_LOG}"

  # Build downstream CLI arguments
  DS_ARGS=( --results-dir "$RESULTS_DIR" --metrics-env "$METRICS_ENV" )

  if [[ "${RUN_MACREL_ATTACH}" -eq 1 ]]; then
    DS_ARGS+=( --macrel-attach-only )
    DS_ARGS+=( --predicted-smorfs "$PREDICTED_SMORFS_TSV" )
    DS_ARGS+=( --id-col "$MACREL_ID_COL" )
    DS_ARGS+=( --seq-col "$MACREL_SEQ_COL" )

    if [[ -n "${MACREL_ATTACH_OUT}" ]]; then
      DS_ARGS+=( --macrel-attach-out "$MACREL_ATTACH_OUT" )
    fi
  fi

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

  if [[ "${RUN_MACREL_ATTACH}" -eq 1 ]]; then
    if [[ -z "${PREDICTED_SMORFS_TSV}" ]]; then
      echo "ERROR: --predicted-smorfs must be provided with --macrel-attach-only"
      exit 1
    fi
    if [[ -z "${MACREL_SEQ_COL}" ]]; then
      echo "ERROR: --seq-col must be provided with --macrel-attach-only"
      exit 1
    fi
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
