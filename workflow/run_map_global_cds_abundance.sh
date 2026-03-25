#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# Script: workflow/run_map_global_cds_abundance.sh
#
# Purpose:
#   1) Build representative GLOBAL CDS FASTA from:
#      results/smorfs/GLOBAL/mmseqs/cluster_map.tsv
#   2) Submit Slurm array jobs to map each library separately
#      back to that representative FASTA
#
# Inputs:
#   - results/smorfs/GLOBAL/mmseqs/cluster_map.tsv
#   - results/smorfs/<ENV>_GLOBAL/catalog/cds_all.fna
#   - metadata/metagenome_files.txt (SampleID<TAB>ABS_FASTQ)
#
# Outputs:
#   results/abundance/global_cds/reference/global_rep_cds.fna
#   results/abundance/global_cds/reference/global_rep_cds.metadata.tsv
#   results/abundance/global_cds/manifests/library_manifest.tsv
#   results/abundance/global_cds/bam/
# ==========================================================

RESULTS_DIR="${RESULTS_DIR:-results}"
METADATA_MAP="${METADATA_MAP:-metadata/metagenome_files.txt}"

PARTITION="${PARTITION:-short}"
TIME="${TIME:-04:00:00}"
CPUS="${CPUS:-8}"
MEM="${MEM:-32G}"

THREADS="${THREADS:-${CPUS}}"

ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ABUND_ENV_NAME="${ABUND_ENV_NAME:-smorf_abundance_env}"

BUILD_REF=1
RUN_MAPPING=1
SAMPLE_ID=""

CLUSTER_MAP="${RESULTS_DIR}/smorfs/GLOBAL/mmseqs/cluster_map.tsv"
OUT_ROOT="${RESULTS_DIR}/abundance/global_cds"
REF_DIR="${OUT_ROOT}/reference"
MANIFEST_DIR="${OUT_ROOT}/manifests"
LOG_DIR="logs"

GLOBAL_REF_FASTA="${REF_DIR}/global_rep_cds.fna"
GLOBAL_REF_META="${REF_DIR}/global_rep_cds.metadata.tsv"
LIB_MANIFEST="${MANIFEST_DIR}/library_manifest.tsv"

die(){ echo "ERROR: $*" >&2; exit 1; }
msg(){ echo "[$(date +'%F %T')] $*" >&2; }

usage() {
  cat <<EOF
Usage:
  bash workflow/run_map_global_cds_abundance.sh [options]

Options:
  --results-dir DIR        (default: results)
  --metadata-map TSV       (default: metadata/metagenome_files.txt)
  --partition NAME         (default: short)
  --time HH:MM:SS          (default: 04:00:00)
  --cpus INT               (default: 8)
  --mem STRING             (default: 32G)
  --sample-id ID           Only map FASTQs belonging to one SampleID
  --build-ref-only         Only build representative CDS FASTA
  --map-only               Skip reference build, only submit mapping jobs
  --env-name NAME          Conda env name/prefix suffix (default: smorf_abundance_env)
  -h, --help               Show this help

Notes:
  - cluster_map.tsv is expected at:
      results/smorfs/GLOBAL/mmseqs/cluster_map.tsv
  - representative IDs are taken as unique values from column 2 of cluster_map.tsv
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --metadata-map) METADATA_MAP="$2"; shift 2 ;;
    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; THREADS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --sample-id) SAMPLE_ID="$2"; shift 2 ;;
    --build-ref-only) BUILD_REF=1; RUN_MAPPING=0; shift ;;
    --map-only) BUILD_REF=0; RUN_MAPPING=1; shift ;;
    --env-name) ABUND_ENV_NAME="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown arg: $1" ;;
  esac
done

CLUSTER_MAP="${RESULTS_DIR}/smorfs/GLOBAL/mmseqs/cluster_map.tsv"
OUT_ROOT="${RESULTS_DIR}/abundance/global_cds"
REF_DIR="${OUT_ROOT}/reference"
MANIFEST_DIR="${OUT_ROOT}/manifests"

GLOBAL_REF_FASTA="${REF_DIR}/global_rep_cds.fna"
GLOBAL_REF_META="${REF_DIR}/global_rep_cds.metadata.tsv"
LIB_MANIFEST="${MANIFEST_DIR}/library_manifest.tsv"

mkdir -p "${REF_DIR}" "${MANIFEST_DIR}" "${LOG_DIR}"

build_reference() {
  [[ -f "${CLUSTER_MAP}" ]] || die "Missing cluster map: ${CLUSTER_MAP}"

  msg "Building GLOBAL representative CDS FASTA via srun"
  srun \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${PWD}" \
    python workflow/build_global_rep_cds_from_cluster_map.py \
      --cluster-map "${CLUSTER_MAP}" \
      --results-dir "${RESULTS_DIR}" \
      --out-fasta "${GLOBAL_REF_FASTA}" \
      --out-meta "${GLOBAL_REF_META}"

  [[ -s "${GLOBAL_REF_FASTA}" ]] || die "Representative CDS FASTA was not created: ${GLOBAL_REF_FASTA}"
  [[ -s "${GLOBAL_REF_META}" ]] || die "Representative CDS metadata TSV was not created: ${GLOBAL_REF_META}"

  msg "Reference FASTA: ${GLOBAL_REF_FASTA}"
  msg "Reference meta : ${GLOBAL_REF_META}"
}

build_manifest() {
  [[ -f "${METADATA_MAP}" ]] || die "Missing metadata map: ${METADATA_MAP}"

  if [[ -n "${SAMPLE_ID}" ]]; then
    awk -F $'\t' -v sid="${SAMPLE_ID}" 'NF>=2 && $1==sid {print $1 "\t" $2}' "${METADATA_MAP}" > "${LIB_MANIFEST}"
  else
    awk -F $'\t' 'NF>=2 {print $1 "\t" $2}' "${METADATA_MAP}" > "${LIB_MANIFEST}"
  fi

  [[ -s "${LIB_MANIFEST}" ]] || die "Manifest is empty: ${LIB_MANIFEST}"
  msg "Library manifest: ${LIB_MANIFEST}"
  msg "Number of libraries: $(wc -l < "${LIB_MANIFEST}")"
}

submit_mapping_array() {
  [[ -s "${GLOBAL_REF_FASTA}" ]] || die "Representative CDS FASTA missing/empty: ${GLOBAL_REF_FASTA}"

  local njobs
  njobs="$(wc -l < "${LIB_MANIFEST}")"
  [[ "${njobs}" -ge 1 ]] || die "No jobs to submit."

  msg "Submitting Slurm array: 1-${njobs}"
  sbatch \
    --job-name=map_global_cds \
    --partition="${PARTITION}" \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --array="1-${njobs}" \
    --output="${LOG_DIR}/map_global_cds_%A_%a.out" \
    --error="${LOG_DIR}/map_global_cds_%A_%a.err" \
    --export=ALL,\
RESULTS_DIR="${RESULTS_DIR}",\
THREADS="${THREADS}",\
ENV_PREFIX_DIR="${ENV_PREFIX_DIR}",\
ABUND_ENV_NAME="${ABUND_ENV_NAME}",\
ABUND_LIBRARY_MANIFEST="${LIB_MANIFEST}",\
ABUND_GLOBAL_REF="${GLOBAL_REF_FASTA}" \
    workflow/map_global_cds_array_task.sh
}

main() {
  msg "============================================"
  msg "GLOBAL CDS abundance mapping runner"
  msg "RESULTS_DIR    : ${RESULTS_DIR}"
  msg "METADATA_MAP   : ${METADATA_MAP}"
  msg "CLUSTER_MAP    : ${CLUSTER_MAP}"
  msg "BUILD_REF      : ${BUILD_REF}"
  msg "RUN_MAPPING    : ${RUN_MAPPING}"
  msg "SAMPLE_ID      : ${SAMPLE_ID:-<all>}"
  msg "PARTITION/TIME : ${PARTITION} / ${TIME}"
  msg "CPUS/MEM       : ${CPUS} / ${MEM}"
  msg "ABUND_ENV_NAME : ${ABUND_ENV_NAME}"
  msg "============================================"

  if [[ "${BUILD_REF}" -eq 1 ]]; then
    build_reference
  fi

  if [[ "${RUN_MAPPING}" -eq 1 ]]; then
    build_manifest
    submit_mapping_array
  fi

  msg "DONE."
}

main "$@" 
