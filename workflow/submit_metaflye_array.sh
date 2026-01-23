#!/usr/bin/env bash
# ==========================================================
# Script: workflow/submit_metaflye_array.sh
# Purpose: Build FASTQ list from data/ and submit Slurm array
#          (1 task per FASTQ) to run MetaFlye in parallel.
# ==========================================================
set -euo pipefail

PARTITION="${PARTITION:-short}"
TIME="${TIME:-12:00:00}"
CPUS="${CPUS:-16}"
MEM="${MEM:-64G}"
WDIR="${WDIR:-$PWD}"

FASTQ_DIR="${FASTQ_DIR:-data}"
RESULTS_DIR="${RESULTS_DIR:-results}"

# Optional knobs exported to tasks
READ_MODE="${READ_MODE:-nano-raw}"     # nano-raw recommended
MIN_OVERLAP="${MIN_OVERLAP:-}"         # empty => auto
GENOME_SIZE="${GENOME_SIZE:-}"         # usually empty for metagenomes

mkdir -p logs metadata

TS="$(date +%Y%m%d_%H%M%S)"
LIST_FILE="metadata/metaflye_fastqs_${TS}.list"

# Build FASTQ list (absolute paths)
shopt -s nullglob
mapfile -t fastqs < <(ls -1 "${FASTQ_DIR}"/*.fastq "${FASTQ_DIR}"/*.fastq.gz "${FASTQ_DIR}"/*.fq "${FASTQ_DIR}"/*.fq.gz 2>/dev/null | sort -V)

if [[ "${#fastqs[@]}" -eq 0 ]]; then
  echo "ERROR: No FASTQ files found in ${FASTQ_DIR}/" >&2
  exit 1
fi

# Write absolute paths
: > "${LIST_FILE}"
for f in "${fastqs[@]}"; do
  # make absolute
  if [[ "${f}" = /* ]]; then
    echo "${f}" >> "${LIST_FILE}"
  else
    echo "${WDIR}/${f}" >> "${LIST_FILE}"
  fi
done

N="${#fastqs[@]}"
echo "MetaFlye: Found ${N} FASTQ file(s)."
echo "FASTQ list written to: ${LIST_FILE}"

OUT_LOG="logs/metaflye_array_${TS}.%A_%a.out"
ERR_LOG="logs/metaflye_array_${TS}.%A_%a.err"

# Submit array: 1..N
jid="$(
sbatch \
  --partition="${PARTITION}" \
  --time="${TIME}" \
  --cpus-per-task="${CPUS}" \
  --mem="${MEM}" \
  --chdir="${WDIR}" \
  --job-name="metaflye" \
  --array="1-${N}" \
  --output="${OUT_LOG}" \
  --error="${ERR_LOG}" \
  --export=ALL,THREADS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",METAFlyE_FASTQ_LIST="${LIST_FILE}",READ_MODE="${READ_MODE}",MIN_OVERLAP="${MIN_OVERLAP}",GENOME_SIZE="${GENOME_SIZE}" \
  workflow/metaflye_array_task.sh \
  | awk '{print $4}'
)"

echo "Submitted MetaFlye array job: ${jid}"
echo "Per-task logs:"
echo "  ${OUT_LOG}"
echo "  ${ERR_LOG}"
