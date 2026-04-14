#!/usr/bin/env bash
set -euo pipefail

PARTITION="${PARTITION:-short}"
TIME="${TIME:-04:00:00}"
CPUS="${CPUS:-8}"
MEM="${MEM:-32G}"
WDIR="${WDIR:-$PWD}"

ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ABUND_ENV_NAME="${ABUND_ENV_NAME:-smorf_aldex2_env}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ABUND_ENV_NAME}"

READ_MAPPING_STATS_ROOT="${READ_MAPPING_STATS_ROOT:-/scratch/t.sousa/data_used/read_mapping/stats}"
OUTDIR="${OUTDIR:-results/differential_abundance/aldex2}"
ENVS="${ENVS:-campina,peneira,uniforme,riparia}"

RUN_FLAGSTAT=0
RUN_PREPARE=0
RUN_ALDEX2=0
CHECK_INSTALL_ONLY=0

MC_SAMPLES="${MC_SAMPLES:-128}"
DENOM="${DENOM:-all}"
MIN_COUNT="${MIN_COUNT:-15}"
MIN_SAMPLES="${MIN_SAMPLES:-5}"
GROUP_COL="${GROUP_COL:-environment}"
USE_MC="${USE_MC:-FALSE}"

COUNTS_TSV=""
METADATA_TSV=""

log(){ echo "[$(date +'%F %T')] $*" >&2; }
die(){ echo "ERROR: $*" >&2; exit 1; }
have_cmd(){ command -v "$1" >/dev/null 2>&1; }

usage() {
  cat <<EOF
Usage:
  bash workflow/run_differential_abundance_aldex2.sh [options]

Modes:
  --flagstat-only
  --prepare-only
  --aldex2-only
  --all
  --check-install-only

Options:
  --partition STR
  --time HH:MM:SS
  --cpus INT
  --mem STR
  --wd PATH
  --env-name STR
  --env-prefix-dir PATH
  --stats-root PATH
  --outdir PATH
  --envs STR
  --counts-tsv PATH
  --metadata-tsv PATH
  --mc-samples INT
  --denom STR
  --min-count INT
  --min-samples INT
  --group-col STR
  --use-mc TRUE|FALSE
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --flagstat-only) RUN_FLAGSTAT=1; shift ;;
    --prepare-only) RUN_PREPARE=1; shift ;;
    --aldex2-only) RUN_ALDEX2=1; shift ;;
    --all) RUN_FLAGSTAT=1; RUN_PREPARE=1; RUN_ALDEX2=1; shift ;;
    --check-install-only) CHECK_INSTALL_ONLY=1; shift ;;

    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --wd) WDIR="$2"; shift 2 ;;

    --env-name)
      ABUND_ENV_NAME="$2"
      ENV_PREFIX="${ENV_PREFIX_DIR}/${ABUND_ENV_NAME}"
      shift 2
      ;;
    --env-prefix-dir)
      ENV_PREFIX_DIR="$2"
      ENV_PREFIX="${ENV_PREFIX_DIR}/${ABUND_ENV_NAME}"
      shift 2
      ;;
    --stats-root) READ_MAPPING_STATS_ROOT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --envs) ENVS="$2"; shift 2 ;;
    --counts-tsv) COUNTS_TSV="$2"; shift 2 ;;
    --metadata-tsv) METADATA_TSV="$2"; shift 2 ;;
    --mc-samples) MC_SAMPLES="$2"; shift 2 ;;
    --denom) DENOM="$2"; shift 2 ;;
    --min-count) MIN_COUNT="$2"; shift 2 ;;
    --min-samples) MIN_SAMPLES="$2"; shift 2 ;;
    --group-col) GROUP_COL="$2"; shift 2 ;;
    --use-mc) USE_MC="$2"; shift 2 ;;

    -h|--help) usage; exit 0 ;;
    *) die "Unknown arg: $1" ;;
  esac
done

if [[ "$RUN_FLAGSTAT" -eq 0 && "$RUN_PREPARE" -eq 0 && "$RUN_ALDEX2" -eq 0 && "$CHECK_INSTALL_ONLY" -eq 0 ]]; then
  RUN_FLAGSTAT=1
  RUN_PREPARE=1
  RUN_ALDEX2=1
fi

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
  die "conda not found"
}

ensure_env_dir() {
  init_conda
  mkdir -p "$ENV_PREFIX_DIR"
}

create_aldex2_env_if_missing() {
  ensure_env_dir

  if [[ -d "$ENV_PREFIX" && -x "$ENV_PREFIX/bin/Rscript" ]]; then
    log "ALDEx2 env already exists: $ENV_PREFIX"
    return 0
  fi

  local installer="conda"
  if have_cmd mamba; then
    installer="mamba"
  fi

  log "Creating dedicated ALDEx2 env at: $ENV_PREFIX"
  "${installer}" create -y -p "$ENV_PREFIX" -c conda-forge -c bioconda \
    r-base r-data.table r-ggplot2 bioconductor-aldex2

  [[ -x "$ENV_PREFIX/bin/Rscript" ]] || die "Failed to create ALDEx2 env at $ENV_PREFIX"
}

activate_env() {
  init_conda
  [[ -d "$ENV_PREFIX" ]] || die "Env prefix not found: $ENV_PREFIX"
  conda activate "$ENV_PREFIX"
}

ensure_aldex2_stack() {
  create_aldex2_env_if_missing
  activate_env

  log "Checking required R packages inside ${ENV_PREFIX}"
  Rscript - <<'RS'
needed <- c("data.table", "ggplot2", "ALDEx2", "edgeR")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  cat("Missing R packages:", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
}
cat("R package stack OK\n")
cat("ALDEx2 version:", as.character(utils::packageVersion("ALDEx2")), "\n")
cat("edgeR version:", as.character(utils::packageVersion("edgeR")), "\n")
RS
}

check_workflow_script() {
  [[ -f "${WDIR}/workflow/aldex2_global_da.R" ]] || die "Missing R workflow script: ${WDIR}/workflow/aldex2_global_da.R"
}

run_r_step() {
  local step="$1"; shift
  check_workflow_script
  activate_env

  srun \
    --partition="$PARTITION" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$CPUS" \
    --mem="$MEM" \
    --time="$TIME" \
    --chdir="$WDIR" \
    Rscript workflow/aldex2_global_da.R --step "$step" "$@"
}

main() {
  log "============================================"
  log "ALDEx2 differential abundance runner"
  log "ENV_PREFIX              : $ENV_PREFIX"
  log "READ_MAPPING_STATS_ROOT : $READ_MAPPING_STATS_ROOT"
  log "OUTDIR                  : $OUTDIR"
  log "ENVS                    : $ENVS"
  log "RUN_FLAGSTAT            : $RUN_FLAGSTAT"
  log "RUN_PREPARE             : $RUN_PREPARE"
  log "RUN_ALDEX2              : $RUN_ALDEX2"
  log "CHECK_INSTALL_ONLY      : $CHECK_INSTALL_ONLY"
  log "MC_SAMPLES              : $MC_SAMPLES"
  log "DENOM                   : $DENOM"
  log "MIN_COUNT               : $MIN_COUNT"
  log "MIN_SAMPLES             : $MIN_SAMPLES"
  log "GROUP_COL               : $GROUP_COL"
  log "USE_MC                  : $USE_MC"
  log "============================================"

  ensure_aldex2_stack

  if [[ "$CHECK_INSTALL_ONLY" -eq 1 ]]; then
    log "Check/install only requested. Exiting."
    exit 0
  fi

  mkdir -p "$OUTDIR"

  if [[ "$RUN_FLAGSTAT" -eq 1 ]]; then
    run_r_step flagstat \
      --stats-root "$READ_MAPPING_STATS_ROOT" \
      --outdir "$OUTDIR" \
      --envs "$ENVS"
  fi

  if [[ "$RUN_PREPARE" -eq 1 ]]; then
    run_r_step prepare \
      --stats-root "$READ_MAPPING_STATS_ROOT" \
      --outdir "$OUTDIR" \
      --envs "$ENVS"
  fi

  if [[ "$RUN_ALDEX2" -eq 1 ]]; then
    [[ -n "$COUNTS_TSV" ]] || COUNTS_TSV="$OUTDIR/aldex2_counts_matrix.tsv"
    [[ -n "$METADATA_TSV" ]] || METADATA_TSV="$OUTDIR/aldex2_sample_metadata.tsv"

    run_r_step run \
      --outdir "$OUTDIR" \
      --counts-tsv "$COUNTS_TSV" \
      --metadata-tsv "$METADATA_TSV" \
      --group-col "$GROUP_COL" \
      --mc-samples "$MC_SAMPLES" \
      --denom "$DENOM" \
      --min-count "$MIN_COUNT" \
      --min-samples "$MIN_SAMPLES" \
      --use-mc "$USE_MC"
  fi

  log "DONE."
}

main "$@"
