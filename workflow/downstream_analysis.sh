#!/usr/bin/env bash
# ==========================================================
# Script: workflow/downstream_analysis.sh
# Purpose:
#   - Create (once) a small conda env to parse MetaFlye flye.log files
#   - Step A (Python): parse flye.log -> metrics TSV
#   - Step B (R): read TSV -> boxplots for each column
#
# Default inputs:
#   results/assembly_metaflye/*/flye.log
#
# Default outputs:
#   results/assembly_metaflye/finalize_metrics.tsv
#   results/assembly_metaflye/finalize_boxplots/boxplot_*.png
# ==========================================================
set -euo pipefail

# ---------------------------
# Defaults
# ---------------------------
ENV_NAME="${METRICS_ENV:-metaflye_metrics_env}"
ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"

BASE_RESULTS_DIR="${RESULTS_DIR:-results}"

OUT_TSV="${OUT_TSV:-${BASE_RESULTS_DIR}/assembly_metaflye/finalize_metrics.tsv}"
PLOTS_DIR="${PLOTS_DIR:-${BASE_RESULTS_DIR}/assembly_metaflye/finalize_boxplots}"

RUN_PARSE=1
RUN_PLOTS=1

# ---------------------------
# Helpers
# ---------------------------
die() { echo "ERROR: $*" >&2; exit 1; }
msg() { echo "[`date +'%F %T'`] $*" >&2; }
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
  mkdir -p "${ENV_PREFIX_DIR}"

  local lockdir="${ENV_PREFIX}.lockdir"

  if [[ -d "${ENV_PREFIX}" ]]; then
    msg "Env exists: ${ENV_PREFIX}"
  else
    while ! mkdir "${lockdir}" 2>/dev/null; do
      msg "Waiting for env lock..."
      sleep 5
    done

    cleanup_lock() { rmdir "${lockdir}" 2>/dev/null || true; }
    trap cleanup_lock EXIT INT TERM

    if [[ -d "${ENV_PREFIX}" ]]; then
      msg "Env appeared while waiting. Continuing."
      return 0
    fi

    msg "Creating conda env at: ${ENV_PREFIX}"
    local installer="conda"
    if have_cmd mamba; then installer="mamba"; fi
    "${installer}" config set channel_priority strict >/dev/null 2>&1 || true

    # Python for parsing + R for plotting
    "${installer}" create -y -p "${ENV_PREFIX}" --override-channels \
      -c conda-forge \
      python=3.10 \
      pandas \
      matplotlib \
      r-base=4.3 \
      r-ggplot2 \
      r-readr \
      r-dplyr

    msg "Env created: ${ENV_PREFIX}"

    trap - EXIT INT TERM
    cleanup_lock
  fi

  set +u
  conda activate "${ENV_PREFIX}"
  set -u

  msg "Python: $(python -V)"
  msg "R:      $(R --version | head -n 1)"
}

step_a_parse() {
  msg "Step A: parse flye.log -> TSV"
  python workflow/summarize_flye_logs.py "${BASE_RESULTS_DIR}" --out "${OUT_TSV}" --no-plots
  msg "Wrote: ${OUT_TSV}"
}

step_b_plots_r() {
  msg "Step B: boxplots (R)"
  Rscript workflow/plot_metaflye_metrics.R "${OUT_TSV}" "${PLOTS_DIR}"
  msg "Plots: ${PLOTS_DIR}/boxplot_*.png"
}

usage() {
  cat <<EOF
Usage:
  bash workflow/downstream_analysis.sh [options]

Options:
  --results-dir PATH     Base results dir (default: results)
  --out PATH             Output TSV (default: results/assembly_metaflye/finalize_metrics.tsv)
  --plots-dir PATH       Output plots dir (default: results/assembly_metaflye/finalize_boxplots)
  --metrics-env STR      Env name under envs/ (default: metaflye_metrics_env)
  --parse-only           Run only Step A (TSV)
  --plots-only           Run only Step B (plots; expects TSV exists)
  -h, --help             Show help

Behavior:
  - On first run: creates the conda env, activates it, then runs A and B.
  - On later runs: reuses the env, then runs A and B (unless you select only one step).
EOF
}

# ---------------------------
# CLI
# ---------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) BASE_RESULTS_DIR="$2"; shift 2 ;;
    --out) OUT_TSV="$2"; shift 2 ;;
    --plots-dir) PLOTS_DIR="$2"; shift 2 ;;
    --metrics-env) ENV_NAME="$2"; ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"; shift 2 ;;
    --parse-only) RUN_PLOTS=0; shift ;;
    --plots-only) RUN_PARSE=0; shift ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1 (use --help)" ;;
  esac
done

# ---------------------------
# Main
# ---------------------------
ensure_env_once

[[ -f workflow/summarize_flye_logs.py ]] || die "Missing: workflow/summarize_flye_logs.py"
[[ -f workflow/plot_metaflye_metrics.R ]] || die "Missing: workflow/plot_metaflye_metrics.R"

if [[ "${RUN_PARSE}" -eq 1 ]]; then
  step_a_parse
fi

if [[ "${RUN_PLOTS}" -eq 1 ]]; then
  step_b_plots_r
fi

msg "DONE"
