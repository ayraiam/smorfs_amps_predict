#!/usr/bin/env bash
set -euo pipefail

PARTITION="short"
TIME="12:00:00"
CPUS="8"
MEM="32G"
WDIR="$PWD"

ENV_NAME="global_cds_annot_env"
INPUT_FNA="results/abundance/global_cds/reference/global_rep_cds.fna"
OUTDIR="results/annotation/global_cds"

SWISSPROT_DB=""
PFAM_DB=""
EGGNOG_DATA_DIR=""

RUN_INIT=0
RUN_ORTHOLOGY=0
RUN_DOMAINS=0
RUN_FUNCTIONAL=0
RUN_TAXONOMY=0
CREATE_ENV=0

usage() {
  cat <<EOF
Usage: bash workflow/run_annotate_global_cds.sh [options]

Core:
  --partition STR
  --time HH:MM:SS
  --cpus INT
  --mem STR
  --wd PATH
  --env-name STR
  --input-fna PATH
  --outdir PATH

Databases:
  --swissprot-db PATH       Pre-built DIAMOND protein DB for orthology/best-hit search
  --pfam-db PATH            Pfam-A.hmm (already hmmpress'ed)
  --eggnog-data-dir PATH    eggNOG-mapper data directory

Modes:
  --create-env
  --all
  --init-only
  --orthology-only
  --domains-only
  --functional-only
  --taxonomy-only

Notes:
  - init creates the translated proteins FASTA and initializes the master TSV.
  - orthology runs DIAMOND blastp against --swissprot-db and updates best-hit columns.
  - domains runs hmmscan against --pfam-db and updates conserved_domains.
  - functional runs eggNOG-mapper and updates enzyme / GO / KEGG / COG columns.
  - taxonomy fills taxonomic_hint, prioritizing DIAMOND species names and then eggNOG tax scope.
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --partition) PARTITION="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --cpus) CPUS="$2"; shift 2 ;;
    --mem) MEM="$2"; shift 2 ;;
    --wd) WDIR="$2"; shift 2 ;;
    --env-name) ENV_NAME="$2"; shift 2 ;;
    --input-fna) INPUT_FNA="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --swissprot-db) SWISSPROT_DB="$2"; shift 2 ;;
    --pfam-db) PFAM_DB="$2"; shift 2 ;;
    --eggnog-data-dir) EGGNOG_DATA_DIR="$2"; shift 2 ;;
    --create-env) CREATE_ENV=1; shift 1 ;;
    --all)
      RUN_INIT=1; RUN_ORTHOLOGY=1; RUN_DOMAINS=1; RUN_FUNCTIONAL=1; RUN_TAXONOMY=1
      shift 1
      ;;
    --init-only) RUN_INIT=1; shift 1 ;;
    --orthology-only) RUN_ORTHOLOGY=1; shift 1 ;;
    --domains-only) RUN_DOMAINS=1; shift 1 ;;
    --functional-only) RUN_FUNCTIONAL=1; shift 1 ;;
    --taxonomy-only) RUN_TAXONOMY=1; shift 1 ;;
    -h|--help) usage ;;
    *) echo "Unknown argument: $1"; usage ;;
  esac
done

mkdir -p "${OUTDIR}"/{01_init,02_orthology,03_domains,04_functional,05_taxonomy,logs}

TABLE_TSV="${OUTDIR}/global_rep_cds_annotation.tsv"
PROTEINS_FAA="${OUTDIR}/01_init/global_rep_cds.faa"
ORTHO_TSV="${OUTDIR}/02_orthology/diamond_best_hits.tsv"
DOMAIN_TBLOUT="${OUTDIR}/03_domains/pfam.domtblout"
EGGNOG_PREFIX="${OUTDIR}/04_functional/global_rep_cds"

if [[ "${CREATE_ENV}" -eq 1 ]]; then
  if command -v mamba >/dev/null 2>&1; then
    PKG_MGR="mamba"
  else
    PKG_MGR="conda"
  fi

  "${PKG_MGR}" create -y -n "${ENV_NAME}" \
    -c conda-forge -c bioconda \
    python=3.11 biopython pandas seqkit diamond hmmer eggnog-mapper
  echo "Created env: ${ENV_NAME}"
  exit 0
fi

run_step() {
  local step="$1"
  local out_log="${OUTDIR}/logs/${step}.out"
  local err_log="${OUTDIR}/logs/${step}.err"

  srun \
    --partition="${PARTITION}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --mem="${MEM}" \
    --time="${TIME}" \
    --chdir="${WDIR}" \
    --export=ALL,ANNOT_ENV_NAME="${ENV_NAME}",INPUT_FNA="${INPUT_FNA}",OUTDIR="${OUTDIR}",TABLE_TSV="${TABLE_TSV}",PROTEINS_FAA="${PROTEINS_FAA}",ORTHO_TSV="${ORTHO_TSV}",DOMAIN_TBLOUT="${DOMAIN_TBLOUT}",EGGNOG_PREFIX="${EGGNOG_PREFIX}",SWISSPROT_DB="${SWISSPROT_DB}",PFAM_DB="${PFAM_DB}",EGGNOG_DATA_DIR="${EGGNOG_DATA_DIR}",THREADS="${CPUS}" \
    /bin/bash workflow/annot_global_cds_job.sh --step "${step}" \
    >>"${out_log}" 2>>"${err_log}"
}

if [[ "${RUN_INIT}" -eq 1 ]]; then
  run_step init
fi
if [[ "${RUN_ORTHOLOGY}" -eq 1 ]]; then
  [[ -n "${SWISSPROT_DB}" ]] || { echo "ERROR: --swissprot-db is required for orthology."; exit 1; }
  run_step orthology
fi
if [[ "${RUN_DOMAINS}" -eq 1 ]]; then
  [[ -n "${PFAM_DB}" ]] || { echo "ERROR: --pfam-db is required for domains."; exit 1; }
  run_step domains
fi
if [[ "${RUN_FUNCTIONAL}" -eq 1 ]]; then
  [[ -n "${EGGNOG_DATA_DIR}" ]] || { echo "ERROR: --eggnog-data-dir is required for functional annotation."; exit 1; }
  run_step functional
fi
if [[ "${RUN_TAXONOMY}" -eq 1 ]]; then
  run_step taxonomy
fi

echo "Done. Master table:"
echo "  ${TABLE_TSV}"
echo "Logs:"
echo "  ${OUTDIR}/logs/"
