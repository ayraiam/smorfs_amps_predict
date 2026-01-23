#!/usr/bin/env bash
# ==========================================================
# Script: workflow/metaflye_array_task.sh
# Purpose: Slurm array task runner: assemble ONE FASTQ using Flye --meta
# Requires:
#   - METAFlyE_FASTQ_LIST env var (path to list of FASTQs; 1 per line)
#   - SLURM_ARRAY_TASK_ID is set
# Output:
#   ${RESULTS_DIR}/assembly_metaflye/<SAMPLE>/
# Logs:
#   handled by sbatch --output/--error in submit script
# ==========================================================
set -euo pipefail

RESULTS_DIR="${RESULTS_DIR:-results}"
THREADS="${THREADS:-8}"

ENV_NAME="${METAFlyE_ENV_NAME:-metaflye}"
ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"

ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-${RESULTS_DIR}/assembly_metaflye}"

# Read type: raw ONT by default (recommended unless we truly have corrected reads)
READ_MODE="${READ_MODE:-nano-raw}"  # nano-raw | nano-hq | nano-corr (use nano-corr ONLY if externally corrected)

# Optional: set MIN_OVERLAP only if we *really* want to override Flye auto-choice
MIN_OVERLAP="${MIN_OVERLAP:-}"      # empty => do NOT pass --min-overlap

# Optional: usually leave empty for metagenomes
GENOME_SIZE="${GENOME_SIZE:-}"

log() { echo "[$(date --iso-8601=seconds)] $*"; }
die() { echo "ERROR: $*" >&2; exit 1; }
have_cmd() { command -v "$1" >/dev/null 2>&1; }

init_conda() {
  if have_cmd conda; then
    # shellcheck disable=SC1090
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi
  for guess in "$HOME/miniforge3" "$HOME/miniconda3" "/opt/conda"; do
    if [[ -f "${guess}/etc/profile.d/conda.sh" ]]; then
      # shellcheck disable=SC1090
      source "${guess}/etc/profile.d/conda.sh"
      return 0
    fi
  done
  die "conda not found. Load conda module or add conda to PATH."
}

ensure_env_once() {
  init_conda
  mkdir -p "${ENV_PREFIX_DIR}" "${ASSEMBLY_ROOT}"

  # Lock to avoid race conditions across array tasks
  local lockdir="${ENV_PREFIX}.lockdir"

  if [[ -d "${ENV_PREFIX}" ]]; then
    log "Env exists: ${ENV_PREFIX}"
  else
    # Acquire lock (atomic mkdir)
    while ! mkdir "${lockdir}" 2>/dev/null; do
      log "Waiting for env lock..."
      sleep 5
    done

    # Double-check after acquiring lock
    if [[ -d "${ENV_PREFIX}" ]]; then
      log "Env appeared while waiting. Continuing."
      rmdir "${lockdir}" || true
    else
      log "Creating conda env at: ${ENV_PREFIX}"
      local installer="conda"
      if have_cmd mamba; then installer="mamba"; fi

      "${installer}" create -y -p "${ENV_PREFIX}" -c conda-forge -c bioconda \
        flye minimap2 samtools python

      log "Env created: ${ENV_PREFIX}"
      rmdir "${lockdir}" || true
    fi
  fi

  conda activate "${ENV_PREFIX}"
  have_cmd flye || die "flye not found after activating env."
  log "Flye: $(flye --version 2>/dev/null || echo 'version unavailable')"
}

pick_fastq_from_list() {
  [[ -n "${METAFlyE_FASTQ_LIST:-}" ]] || die "METAFlyE_FASTQ_LIST env var not set"
  [[ -f "${METAFlyE_FASTQ_LIST}" ]] || die "FASTQ list not found: ${METAFlyE_FASTQ_LIST}"
  [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] || die "SLURM_ARRAY_TASK_ID not set (are you running as an array task?)"

  local fq
  fq="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${METAFlyE_FASTQ_LIST}")"
  [[ -n "${fq}" ]] || die "No FASTQ found at line ${SLURM_ARRAY_TASK_ID} of ${METAFlyE_FASTQ_LIST}"
  [[ -f "${fq}" ]] || die "FASTQ path does not exist: ${fq}"
  echo "${fq}"
}

run_flye_one() {
  local fq="$1"
  local base sample outdir

  base="$(basename "$fq")"
  sample="$base"
  sample="${sample%.fastq.gz}"
  sample="${sample%.fq.gz}"
  sample="${sample%.fastq}"
  sample="${sample%.fq}"

  outdir="${ASSEMBLY_ROOT}/${sample}"
  mkdir -p "${outdir}"

  # Skip if already assembled
  if [[ -s "${outdir}/assembly.fasta" || -s "${outdir}/assembly.fna" ]]; then
    log "Existing assembly found for ${sample} in ${outdir}. Skipping."
    return 0
  fi

  local -a cmd
  case "${READ_MODE}" in
    nano-raw) cmd=(flye --meta --nano-raw "${fq}") ;;
    nano-hq)  cmd=(flye --meta --nano-hq  "${fq}") ;;
    nano-corr) cmd=(flye --meta --nano-corr "${fq}") ;;
    *) die "Unknown READ_MODE=${READ_MODE}. Use nano-raw|nano-hq|nano-corr" ;;
  esac

  cmd+=(--out-dir "${outdir}" --threads "${THREADS}")

  if [[ -n "${MIN_OVERLAP}" ]]; then
    cmd+=(--min-overlap "${MIN_OVERLAP}")
  fi
  if [[ -n "${GENOME_SIZE}" ]]; then
    cmd+=(--genome-size "${GENOME_SIZE}")
  fi

  log "------------------------------------------------------------"
  log "Sample : ${sample}"
  log "Input  : ${fq}"
  log "Outdir : ${outdir}"
  log "Threads: ${THREADS}"
  log "Cmd    : ${cmd[*]}"

  "${cmd[@]}"

  log "DONE: ${sample}"
}

main() {
  log "Starting MetaFlye array task"
  log "RESULTS_DIR : ${RESULTS_DIR}"
  log "ASSEMBLY_ROOT: ${ASSEMBLY_ROOT}"
  log "THREADS     : ${THREADS}"
  log "ENV_PREFIX  : ${ENV_PREFIX}"
  log "READ_MODE   : ${READ_MODE}"
  log "MIN_OVERLAP : ${MIN_OVERLAP:-<auto>}"

  ensure_env_once

  fq="$(pick_fastq_from_list)"
  run_flye_one "${fq}"

  log "Task finished."
}

main "$@"
