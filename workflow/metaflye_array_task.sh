#!/usr/bin/env bash
# ==========================================================
# Script: workflow/metaflye_array_task.sh
# Purpose: Slurm array task runner: co-assemble ONE SampleID using Flye --meta
# Requires:
#   - METAFlyE_SAMPLE_LIST env var (SampleIDs; 1 per line)
#   - METAFlyE_SAMPLE_MAP  env var (TSV: SampleID<TAB>ABS_FASTQ)
#   - SLURM_ARRAY_TASK_ID is set
# Output:
#   ${RESULTS_DIR}/assembly_metaflye/<SampleID>/
# ==========================================================
set -euo pipefail

RESULTS_DIR="${RESULTS_DIR:-results}"
THREADS="${THREADS:-8}"

ENV_NAME="${METAFlyE_ENV_NAME:-metaflye}"
ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"

ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-${RESULTS_DIR}/assembly_metaflye}"

# Store the huge merged co-assembly FASTQ on scratch (and delete after Flye)
SCRATCH_READS_ROOT="${SCRATCH_READS_ROOT:-/scratch/t.sousa/data_used/metaflye_tmp}"
CLEANUP_COASSEMBLY_READS="${CLEANUP_COASSEMBLY_READS:-1}"

READ_MODE="${READ_MODE:-nano-raw}"  # nano-raw | nano-hq | nano-corr
MIN_OVERLAP="${MIN_OVERLAP:-}"      # empty => do NOT pass --min-overlap
GENOME_SIZE="${GENOME_SIZE:-}"      # usually empty for metagenomes

log() { echo "[$(date --iso-8601=seconds)] $*"; }
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
  mkdir -p "${ENV_PREFIX_DIR}" "${ASSEMBLY_ROOT}"

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
        flye minimap2 samtools pigz python

      log "Env created: ${ENV_PREFIX}"
      rmdir "${lockdir}" || true
    fi
  fi

  conda activate "${ENV_PREFIX}"
  have_cmd flye || die "flye not found after activating env."
  have_cmd pigz || die "pigz not found after activating env."
  log "Flye: $(flye --version 2>/dev/null || echo 'version unavailable')"
}

pick_sample_from_list() {
  [[ -n "${METAFlyE_SAMPLE_LIST:-}" ]] || die "METAFlyE_SAMPLE_LIST env var not set"
  [[ -f "${METAFlyE_SAMPLE_LIST}" ]] || die "Sample list not found: ${METAFlyE_SAMPLE_LIST}"
  [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "SLURM_ARRAY_TASK_ID not set"

  local sid
  sid="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${METAFlyE_SAMPLE_LIST}")"
  [[ -n "${sid}" ]] || die "No SampleID at line ${SLURM_ARRAY_TASK_ID} of ${METAFlyE_SAMPLE_LIST}"
  echo "${sid}"
}

get_fastqs_for_sample() {
  local sid="$1"
  [[ -n "${METAFlyE_SAMPLE_MAP:-}" ]] || die "METAFlyE_SAMPLE_MAP env var not set"
  [[ -f "${METAFlyE_SAMPLE_MAP}" ]] || die "Sample map not found: ${METAFlyE_SAMPLE_MAP}"

  # Print one FASTQ per line
  awk -v s="$sid" -F $'\t' '$1==s {print $2}' "${METAFlyE_SAMPLE_MAP}" | sed '/^$/d'
}

build_coassembly_reads() {
  # Build a single gzipped FASTQ for this SampleID (once), from multiple library FASTQs.
  local sid="$1"
  local outdir="$2"

  mkdir -p "${SCRATCH_READS_ROOT}/${sid}"
  local reads_gz="${SCRATCH_READS_ROOT}/${sid}/reads.coassembly.fastq.gz"
  local manifest="${SCRATCH_READS_ROOT}/${sid}/inputs.fastq.list"

  if [[ -s "${reads_gz}" ]]; then
    log "Co-assembly reads already exist: ${reads_gz}"
    return 0
  fi

  mapfile -t fqs < <(get_fastqs_for_sample "$sid")
  [[ "${#fqs[@]}" -ge 1 ]] || die "No FASTQs found in map for SampleID=${sid}"

  printf "%s\n" "${fqs[@]}" > "${manifest}"

  log "Building co-assembly FASTQ for SampleID=${sid}"
  log "Input FASTQs: ${#fqs[@]}"
  log "Manifest: ${manifest}"
  log "Output: ${reads_gz}"

  # Helper: return 0 if this FASTQ needs per-file duplicate header fixing
  needs_dup_fix() {
    local base
    base="$(basename "$1")"
    case "$base" in
      nanopore_shotgun_RDS26_L02-2500-low_04.fastq.gz| \
      nanopore_shotgun_RDS26_L02-2500-normal_12.fastq.gz| \
      nanopore_shotgun_RDS26_L02-4500-low_07.fastq.gz| \
      nanopore_shotgun_RDS26_L02-4500-normal_15.fastq.gz| \
      nanopore_shotgun_RDS26_LO2-500-low_03.fastq.gz| \
      nanopore_shotgun_RDS26_LO2-500-normal_11.fastq.gz)
        return 0
        ;;
      *)
        return 1
        ;;
    esac
  }

  # Rewrite duplicated read IDs (within a single file stream) by appending _dupN on duplicates only.
  # Preserves everything after the first token on the header line.
  rewrite_dup_headers_awk='
    NR%4==1 {
      line = $0
      first = $1                         # first token, includes leading "@"
      id = first
      sub(/^@/, "", id)                  # strip "@"
      c[id]++
      new = id
      if (c[id] > 1) new = id "_dup" c[id]
      rest = substr(line, length(first)+1)  # keeps leading space if present, or empty
      print "@" new rest
      next
    }
    { print }
  '

  local tmp="${reads_gz}.tmp.$$"

  (
    for f in "${fqs[@]}"; do
      [[ -f "$f" ]] || die "Missing FASTQ: $f"

      if needs_dup_fix "$f"; then
        log "FOUND dup-header file: $(basename "$f") -> rewriting duplicated read IDs with _dupN"
        if [[ "$f" =~ \.gz$ ]]; then
          pigz -dc "$f" | awk "${rewrite_dup_headers_awk}"
        else
          cat "$f" | awk "${rewrite_dup_headers_awk}"
        fi
      else
        if [[ "$f" =~ \.gz$ ]]; then
          pigz -dc "$f"
        else
          cat "$f"
        fi
      fi
    done
  ) | pigz -p "${THREADS}" -c > "${tmp}"

  mv "${tmp}" "${reads_gz}"
  log "Co-assembly FASTQ created: ${reads_gz}"
}

run_flye_sample() {
  local sid="$1"
  local outdir="${ASSEMBLY_ROOT}/${sid}"
  mkdir -p "${outdir}"

  # Skip if already assembled
  if [[ -s "${outdir}/assembly.fasta" || -s "${outdir}/assembly.fna" ]]; then
    log "Existing assembly found for ${sid} in ${outdir}. Skipping."
    return 0
  fi

  build_coassembly_reads "${sid}" "${outdir}"
  local fq="${SCRATCH_READS_ROOT}/${sid}/reads.coassembly.fastq.gz"

  local -a cmd
  case "${READ_MODE}" in
    nano-raw) cmd=(flye --meta --nano-raw "${fq}") ;;
    nano-hq)  cmd=(flye --meta --nano-hq  "${fq}") ;;
    nano-corr) cmd=(flye --meta --nano-corr "${fq}") ;;
    *) die "Unknown READ_MODE=${READ_MODE}. Use nano-raw|nano-hq|nano-corr" ;;
  esac

  cmd+=(--out-dir "${outdir}" --threads "${THREADS}")

  if [[ -n "${MIN_OVERLAP}" ]]; then cmd+=(--min-overlap "${MIN_OVERLAP}"); fi
  if [[ -n "${GENOME_SIZE}" ]]; then cmd+=(--genome-size "${GENOME_SIZE}"); fi

  log "------------------------------------------------------------"
  log "SampleID: ${sid}"
  log "Outdir  : ${outdir}"
  log "Reads   : ${fq}"
  log "Threads : ${THREADS}"
  log "Cmd     : ${cmd[*]}"

  "${cmd[@]}"
  if [[ "${CLEANUP_COASSEMBLY_READS}" -eq 1 ]]; then
    log "Cleaning up scratch co-assembly reads for ${sid}"
    rm -f "${SCRATCH_READS_ROOT}/${sid}/reads.coassembly.fastq.gz" \
          "${SCRATCH_READS_ROOT}/${sid}/inputs.fastq.list" || true
    rmdir "${SCRATCH_READS_ROOT}/${sid}" 2>/dev/null || true
	fi
  log "DONE: ${sid}"
}

main() {
  log "Starting MetaFlye co-assembly array task"
  log "RESULTS_DIR         : ${RESULTS_DIR}"
  log "ASSEMBLY_ROOT       : ${ASSEMBLY_ROOT}"
  log "THREADS             : ${THREADS}"
  log "READ_MODE           : ${READ_MODE}"
  log "MIN_OVERLAP         : ${MIN_OVERLAP:-<auto>}"
  log "METAFlyE_SAMPLE_LIST: ${METAFlyE_SAMPLE_LIST:-<unset>}"
  log "METAFlyE_SAMPLE_MAP : ${METAFlyE_SAMPLE_MAP:-<unset>}"

  ensure_env_once

  sid="$(pick_sample_from_list)"
  run_flye_sample "${sid}"

  log "Task finished."
}

main "$@"
