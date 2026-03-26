#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# Script: workflow/map_global_cds_array_task.sh
# Purpose: Slurm array task runner: map ONE library FASTQ
#          back to the GLOBAL representative CDS FASTA.
#
# Requires env vars:
#   - ABUND_LIBRARY_MANIFEST : TSV with SampleID<TAB>ABS_FASTQ
#   - ABUND_GLOBAL_REF       : representative CDS FASTA
#   - SLURM_ARRAY_TASK_ID
#
# Outputs:
#   /scratch/t.sousa/data_used/read_mapping/bam/<SampleID>/<library>.vs_global_rep_cds.primary_q20.bam
#   /scratch/t.sousa/data_used/read_mapping/bam/<SampleID>/<library>.vs_global_rep_cds.primary_q20.bam.bai
#   /scratch/t.sousa/data_used/read_mapping/stats/<SampleID>/<library>.primary_q20.flagstat.txt
#   /scratch/t.sousa/data_used/read_mapping/stats/<SampleID>/<library>.primary_q20.idxstats.tsv
# ==========================================================

RESULTS_DIR="${RESULTS_DIR:-results}"
THREADS="${THREADS:-8}"

ABUND_ENV_NAME="${ABUND_ENV_NAME:-smorf_abundance_env}"
ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ABUND_ENV_NAME}"

# All mapping outputs go here
READ_MAPPING_ROOT="${READ_MAPPING_ROOT:-/scratch/t.sousa/data_used/read_mapping}"

BAM_ROOT="${READ_MAPPING_ROOT}/bam"
STATS_ROOT="${READ_MAPPING_ROOT}/stats"
TMP_ROOT="${READ_MAPPING_ROOT}/tmp"

log() { echo "[$(date --iso-8601=seconds)] $*" >&2; }
die() { echo "ERROR: $*" >&2; exit 1; }
have_cmd() { command -v "$1" >/dev/null 2>&1; }

init_conda() {
  if have_cmd conda; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi
  for guess in "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "${guess}/etc/profile.d/conda.sh" ]]; then
      source "${guess}/etc/profile.d/conda.sh"
      return 0
    fi
  done
  die "conda not found. Load conda module or add conda to PATH."
}

ensure_env_once() {
  init_conda
  mkdir -p "${ENV_PREFIX_DIR}" "${BAM_ROOT}" "${STATS_ROOT}" "${TMP_ROOT}"

  local lockdir="${ENV_PREFIX}.lockdir"

  if [[ -d "${ENV_PREFIX}" ]]; then
    log "Env exists: ${ENV_PREFIX}"
  else
    while ! mkdir "${lockdir}" 2>/dev/null; do
      log "Waiting for env lock..."
      sleep 5
    done

    if [[ -d "${ENV_PREFIX}" ]]; then
      log "Env appeared while waiting. Continuing."
      rmdir "${lockdir}" || true
    else
      log "Creating conda env at: ${ENV_PREFIX}"
      local installer="conda"
      if have_cmd mamba; then installer="mamba"; fi

      "${installer}" create -y -p "${ENV_PREFIX}" -c conda-forge -c bioconda \
        minimap2 samtools python pandas pysam pigz

      log "Env created: ${ENV_PREFIX}"
      rmdir "${lockdir}" || true
    fi
  fi

  conda activate "${ENV_PREFIX}"
  have_cmd minimap2 || die "minimap2 not found after activating env."
  have_cmd samtools || die "samtools not found after activating env."
  python - <<'PY'
import pandas, pysam
print("[INFO] pandas OK")
print("[INFO] pysam OK")
PY
}

pick_line_from_manifest() {
  [[ -n "${ABUND_LIBRARY_MANIFEST:-}" ]] || die "ABUND_LIBRARY_MANIFEST env var not set"
  [[ -f "${ABUND_LIBRARY_MANIFEST}" ]] || die "Manifest not found: ${ABUND_LIBRARY_MANIFEST}"
  [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "SLURM_ARRAY_TASK_ID not set"

  sed -n "${SLURM_ARRAY_TASK_ID}p" "${ABUND_LIBRARY_MANIFEST}"
}

main() {
  [[ -n "${ABUND_GLOBAL_REF:-}" ]] || die "ABUND_GLOBAL_REF env var not set"
  [[ -f "${ABUND_GLOBAL_REF}" ]] || die "Representative CDS FASTA not found: ${ABUND_GLOBAL_REF}"

  ensure_env_once

  local line sample_id fastq lib_base sample_bam_dir sample_stats_dir filtered_bam tmp_sorted
  line="$(pick_line_from_manifest)"
  [[ -n "${line}" ]] || die "No line found in manifest for array index ${SLURM_ARRAY_TASK_ID}"

  sample_id="$(printf '%s\n' "${line}" | cut -f1)"
  fastq="$(printf '%s\n' "${line}" | cut -f2)"

  [[ -n "${sample_id}" ]] || die "Empty sample_id in manifest line: ${line}"
  [[ -n "${fastq}" ]] || die "Empty fastq in manifest line: ${line}"
  [[ -f "${fastq}" ]] || die "FASTQ not found: ${fastq}"

  lib_base="$(basename "${fastq}")"
  lib_base="${lib_base%.gz}"
  lib_base="${lib_base%.fastq}"
  lib_base="${lib_base%.fq}"

  sample_bam_dir="${BAM_ROOT}/${sample_id}"
  sample_stats_dir="${STATS_ROOT}/${sample_id}"
  mkdir -p "${sample_bam_dir}" "${sample_stats_dir}" "${TMP_ROOT}/${sample_id}"

  filtered_bam="${sample_bam_dir}/${lib_base}.vs_global_rep_cds.primary_q20.bam"
  tmp_sorted="${TMP_ROOT}/${sample_id}/${lib_base}.tmp.sorted.bam"

  if [[ -s "${filtered_bam}" && -s "${filtered_bam}.bai" ]]; then
    log "Filtered BAM already exists for ${sample_id} / ${lib_base}. Skipping."
    exit 0
  fi

  log "------------------------------------------------------------"
  log "SampleID          : ${sample_id}"
  log "FASTQ             : ${fastq}"
  log "LIB_BASE          : ${lib_base}"
  log "REF               : ${ABUND_GLOBAL_REF}"
  log "READ_MAPPING_ROOT : ${READ_MAPPING_ROOT}"
  log "THREADS           : ${THREADS}"

  minimap2 -ax map-ont -t "${THREADS}" "${ABUND_GLOBAL_REF}" "${fastq}" \
    | samtools view -@ "${THREADS}" -bS - \
    | samtools sort -@ "${THREADS}" -o "${tmp_sorted}"

  # Keep primary alignments with MAPQ >= 20
  samtools view -@ "${THREADS}" -b -F 2308 -q 20 "${tmp_sorted}" > "${filtered_bam}"
  samtools index -@ "${THREADS}" "${filtered_bam}"

  samtools flagstat -@ "${THREADS}" "${filtered_bam}" > "${sample_stats_dir}/${lib_base}.primary_q20.flagstat.txt"
  samtools idxstats "${filtered_bam}" > "${sample_stats_dir}/${lib_base}.primary_q20.idxstats.tsv"

  rm -f "${tmp_sorted}"

  log "DONE: ${sample_id} / ${lib_base}"
}

main "$@"
