#!/usr/bin/env bash
set -euo pipefail

ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
REFINE_ENV="${REFINE_EUKS_ENV:-refine_annot_smorf_euks_env}"
CPUS=8
MODE=""
SAMPLE=""
SAMPLES_FILE=""
INPUT_TSV=""
RUN_STEP1=1
RUN_STEP2=1
RUN_STEP3=1

die() { echo "ERROR: $*" >&2; exit 1; }
msg() { echo "[$(date +'%F %T')] $*" >&2; }
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
  die "conda not found."
}

ENV_PREFIX="${ENV_PREFIX_DIR}/${REFINE_ENV}"

create_env() {
  init_conda
  mkdir -p "${ENV_PREFIX_DIR}"

  if [[ -d "${ENV_PREFIX}" ]]; then
    msg "Env exists: ${ENV_PREFIX}"
    return 0
  fi

  local installer="conda"
  have_cmd mamba && installer="mamba"

  "${installer}" create -y -p "${ENV_PREFIX}" --override-channels \
    -c conda-forge -c bioconda \
    python=3.10 pandas biopython gffutils pyranges intervaltree
}

activate_env() {
  [[ -x "${ENV_PREFIX}/bin/python" ]] || die "Env python not found: ${ENV_PREFIX}/bin/python"
  export PATH="${ENV_PREFIX}/bin:${PATH}"
  hash -r
}

run_one_sample() {
  local sample_id="$1"
  local results_dir="${RESULTS_DIR:-results}"
  local sample_dir="${results_dir}/smorfs/${sample_id}"

  [[ -d "${sample_dir}" ]] || die "Sample smorfs dir not found: ${sample_dir}"

  local in_tsv
  if [[ -n "${INPUT_TSV}" ]]; then
    in_tsv="${INPUT_TSV}"
  else
    if [[ -f "${sample_dir}/catalog/predicted_smorfs.with_macrel.refined_euks.tsv" ]]; then
      in_tsv="${sample_dir}/catalog/predicted_smorfs.with_macrel.refined_euks.tsv"
    elif [[ -f "${sample_dir}/catalog/predicted_smorfs.refined_euks.tsv" ]]; then
      in_tsv="${sample_dir}/catalog/predicted_smorfs.refined_euks.tsv"
    elif [[ -f "${sample_dir}/catalog/predicted_smorfs.with_macrel.tsv" ]]; then
      in_tsv="${sample_dir}/catalog/predicted_smorfs.with_macrel.tsv"
    else
      in_tsv="${sample_dir}/catalog/predicted_smorfs.tsv"
    fi
  fi
  [[ -f "${in_tsv}" ]] || die "Input TSV not found: ${in_tsv}"

  local metaeuk_gff="${sample_dir}/fungi/metaeuk/metaeuk_preds.gff"
  [[ -f "${metaeuk_gff}" ]] || die "Missing MetaEuk GFF: ${metaeuk_gff}"

  local fungi_fa="${sample_dir}/contigs/fungi_contigs.fasta"
  [[ -f "${fungi_fa}" ]] || die "Missing fungi contigs FASTA: ${fungi_fa}"

  local cluster_map="${results_dir}/smorfs/GLOBAL/mmseqs/cluster_map.tsv"
  local env_label=""
  if [[ "${sample_id}" =~ ^(RIPARIA|PENEIRA|CAMPINA|UNIFORME)_GLOBAL$ ]]; then
    env_label="${BASH_REMATCH[1]}"
  fi

  local cluster_args=()
  if [[ -f "${cluster_map}" && -n "${env_label}" ]]; then
    cluster_args=( --cluster-map "${cluster_map}" --env-label "${env_label}" )
  fi

  local base="$(basename "${in_tsv}")"
  base="${base%.refined_euks.tsv}"
  base="${base%.tsv}"

  local out_tsv="${sample_dir}/catalog/${base}.refined_euks.tsv"

  python workflow/refine_euks.py \
    --sample "${sample_id}" \
    --input-tsv "${in_tsv}" \
    --metaeuk-gff "${metaeuk_gff}" \
    --fungi-contigs "${fungi_fa}" \
    --out "${out_tsv}" \
    --run-step1 "${RUN_STEP1}" \
    --run-step2 "${RUN_STEP2}" \
    --run-step3 "${RUN_STEP3}" \
    "${cluster_args[@]}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --create-env) MODE="create-env"; shift ;;
    --run) MODE="run"; shift ;;
    --refine-env) REFINE_ENV="$2"; ENV_PREFIX="${ENV_PREFIX_DIR}/${REFINE_ENV}"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --samples-file) SAMPLES_FILE="$2"; shift 2 ;;
    --input-tsv) INPUT_TSV="$2"; shift 2 ;;
    --run-step1) RUN_STEP1="$2"; shift 2 ;;
    --run-step2) RUN_STEP2="$2"; shift 2 ;;
    --run-step3) RUN_STEP3="$2"; shift 2 ;;
    -h|--help) exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -n "${MODE}" ]] || die "Provide --create-env or --run"

if [[ "${MODE}" == "create-env" ]]; then
  create_env
  exit 0
fi

create_env
activate_env

if [[ -n "${SAMPLE}" ]]; then
  run_one_sample "${SAMPLE}"
elif [[ -n "${SAMPLES_FILE}" ]]; then
  while read -r sid; do
    [[ -z "${sid}" ]] && continue
    run_one_sample "${sid}"
  done < "${SAMPLES_FILE}"
else
  die "Provide --sample or --samples-file"
fi
