#!/usr/bin/env bash
set -euo pipefail

PARTITION="${PARTITION:-short}"
TIME="${TIME:-04:00:00}"
CPUS="${CPUS:-8}"
MEM="${MEM:-32G}"
WDIR="${WDIR:-$PWD}"

ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ABUND_ENV_NAME="${ABUND_ENV_NAME:-smorf_abundance_env}"
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
MIN_COUNT="${MIN_COUNT:-10}"
MIN_SAMPLES="${MIN_SAMPLES:-2}"
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
  --envs STR                  Comma-separated, e.g. campina,peneira,uniforme,riparia
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

ensure_env_prefix_exists() {
  [[ -d "$ENV_PREFIX" ]] || die "Env prefix not found: $ENV_PREFIX"
}

conda_install_into_env() {
  local installer="conda"
  if have_cmd mamba; then
    installer="mamba"
  fi

  log "Installing into env: $ENV_PREFIX"
  log "Packages: $*"

  "${installer}" install -y -p "$ENV_PREFIX" -c conda-forge -c bioconda "$@"
}

activate_env() {
  init_conda
  ensure_env_prefix_exists
  conda activate "$ENV_PREFIX"
}

ensure_rscript() {
  activate_env
  if have_cmd Rscript; then
    log "Rscript found in env"
    return 0
  fi

  log "Rscript not found in env. Installing r-base..."
  conda_install_into_env r-base
  conda activate "$ENV_PREFIX"
  have_cmd Rscript || die "Rscript still not found after installing r-base into $ENV_PREFIX"
  log "Rscript installed successfully"
}

ensure_r_packages_cran() {
  ensure_rscript

  log "Checking/installing CRAN R packages needed by aldex2_global_da.R"
  Rscript - <<'RS'
options(repos = c(CRAN = "https://cloud.r-project.org"))

needed <- c("optparse", "data.table", "ggplot2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (length(missing) > 0) {
  install.packages(missing)
}

still_missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing) > 0) {
  cat("Missing CRAN packages after install:", paste(still_missing, collapse=", "), "\n")
  quit(status = 1)
}

cat("CRAN packages OK\n")
RS
}

ensure_aldex2() {
  ensure_r_packages_cran

  log "Checking ALDEx2 inside ${ENV_PREFIX}"
  if Rscript -e 'quit(status = !requireNamespace("ALDEx2", quietly = TRUE))'; then
    log "ALDEx2 already installed"
    Rscript -e 'cat("ALDEx2 version:", as.character(utils::packageVersion("ALDEx2")), "\n")'
    return 0
  fi

  log "ALDEx2 not found. Installing from Bioconductor..."
  Rscript - <<'RS'
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("ALDEx2", ask = FALSE, update = FALSE)

if (!requireNamespace("ALDEx2", quietly = TRUE)) {
  quit(status = 1)
}

cat("ALDEx2 installed successfully\n")
cat("ALDEx2 version:", as.character(utils::packageVersion("ALDEx2")), "\n")
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

  ensure_aldex2

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
