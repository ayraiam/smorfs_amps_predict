#!/usr/bin/env bash
set -euo pipefail

# Defaults
ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
REFINE_ENV="${REFINE_BACS_ENV:-refine_annot_smorf_bacs_env}"
CPUS=8
MODE=""
SAMPLE=""
SAMPLES_FILE=""
INPUT_TSV=""

die() { echo "ERROR: $*" >&2; exit 1; }
msg() { echo "[$(date +'%F %T')] $*" >&2; }
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
  die "conda not found."
}

usage() {
  cat <<EOF
Usage:
  workflow/run_refine_annot_smorf_bacs.sh --create-env --refine-env NAME
  workflow/run_refine_annot_smorf_bacs.sh --run --sample SAMPLE --cpus 8
  workflow/run_refine_annot_smorf_bacs.sh --run --samples-file FILE --cpus 8
Optional:
  --input-tsv PATH    (advanced override; usually leave empty)
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --create-env) MODE="create-env"; shift ;;
    --run) MODE="run"; shift ;;
    --refine-env) REFINE_ENV="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --samples-file) SAMPLES_FILE="$2"; shift 2 ;;
    --input-tsv) INPUT_TSV="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${MODE}" ]] || { usage; exit 1; }

ENV_PREFIX="${ENV_PREFIX_DIR}/${REFINE_ENV}"

create_env() {
  init_conda
  mkdir -p "${ENV_PREFIX_DIR}"

  if [[ -d "${ENV_PREFIX}" ]]; then
    msg "Env exists: ${ENV_PREFIX}"
    return 0
  fi

  msg "Creating refine env at: ${ENV_PREFIX}"
  local installer="conda"
  have_cmd mamba && installer="mamba"

  "${installer}" create -y -p "${ENV_PREFIX}" --override-channels \
    -c conda-forge -c bioconda \
    python=3.10 \
    pandas \
    numpy \
    pyarrow \
    biopython \
    gffutils \
    pyranges \
    intervaltree \
    seqkit

  msg "Env created: ${ENV_PREFIX}"
}

activate_env() {
  init_conda
  # shellcheck disable=SC1090
  conda activate "${ENV_PREFIX}"
  msg "Python: $(python -V)"
}

run_one_sample() {
  local sample_id="$1"
  local results_dir="${RESULTS_DIR:-results}"
  local sample_dir="${results_dir}/smorfs/${sample_id}"
  [[ -d "${sample_dir}" ]] || die "Sample smorfs dir not found: ${sample_dir}"

  # Default input TSV: you can adjust to your actual naming
  local in_tsv
  if [[ -n "${INPUT_TSV}" ]]; then
    in_tsv="${INPUT_TSV}"
  else
    # prefer with_macrel if exists; fallback to predicted_smorfs.tsv
    if [[ -f "${sample_dir}/catalog/predicted_smorfs.with_macrel.tsv" ]]; then
      in_tsv="${sample_dir}/catalog/predicted_smorfs.with_macrel.tsv"
    else
      in_tsv="${sample_dir}/catalog/predicted_smorfs.tsv"
    fi
  fi
  [[ -f "${in_tsv}" ]] || die "Input TSV not found: ${in_tsv}"

  # Prodigal GFF (bacterial)
  local gff="${sample_dir}/bac/prodigal/bac.genes.gff"
  local gff_args=()

  if [[ -f "${gff}" ]]; then
    gff_args=( --prodigal-gff "${gff}" )
  else
    msg "[${sample_id}] WARNING: Prodigal GFF not found; Step2 (embedded) will be skipped: ${gff}"
  fi

  # SmORFinder GFF (needed for smorfinder coordinate mapping + overlap)
  local smorf_gff="${sample_dir}/bac/smorfinder/smorf_output/smorf_output.gff"
  local smorf_args=()

  if [[ -f "${smorf_gff}" ]]; then
    smorf_args=( --smorfinder-gff "${smorf_gff}" )
  else
    msg "[${sample_id}] WARNING: SmORFinder GFF not found; smorfinder mapping/overlap will be skipped: ${smorf_gff}"
  fi

  # Contigs fasta (bacterial)
  local bac_fa="${sample_dir}/contigs/bac_contigs.fasta"
  [[ -f "${bac_fa}" ]] || msg "[${sample_id}] NOTE: bac_contigs.fasta missing: ${bac_fa}"

  local out_tsv="${sample_dir}/catalog/$(basename "${in_tsv%.tsv}").refined_bacs.tsv"

  msg "[${sample_id}] Running refine python"
  python workflow/refine_bacs.py \
    --sample "${sample_id}" \
    --input-tsv "${in_tsv}" \
    "${gff_args[@]}" \
    "${smorf_args[@]}" \
    --bac-contigs "${bac_fa}" \
    --results-dir "${results_dir}" \
    --out "${out_tsv}" \
    --cpus "${CPUS}"
}

if [[ "${MODE}" == "create-env" ]]; then
  create_env
  exit 0
fi

# MODE=run
create_env
activate_env

if [[ -n "${SAMPLE}" ]]; then
  run_one_sample "${SAMPLE}"
elif [[ -n "${SAMPLES_FILE}" ]]; then
  [[ -f "${SAMPLES_FILE}" ]] || die "Samples file not found: ${SAMPLES_FILE}"
  while read -r sid; do
    [[ -z "${sid}" ]] && continue
    run_one_sample "${sid}"
  done < "${SAMPLES_FILE}"
else
  die "Provide --sample or --samples-file"
fi
