#!/usr/bin/env bash
set -euo pipefail

STEP=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --step) STEP="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

[[ -n "${STEP}" ]] || { echo "ERROR: --step is required"; exit 1; }

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found in PATH"
  exit 1
fi

eval "$(conda shell.bash hook)"
conda activate "${ANNOT_ENV_NAME}"

mkdir -p "${OUTDIR}"/{01_init,02_orthology,03_domains,04_functional,05_taxonomy}

case "${STEP}" in
  init)
    python workflow/update_global_annotation_table.py init \
      --input-fna "${INPUT_FNA}" \
      --proteins-faa "${PROTEINS_FAA}" \
      --output-tsv "${TABLE_TSV}"
    ;;

  orthology)
    [[ -s "${PROTEINS_FAA}" ]] || { echo "ERROR: proteins FASTA missing: ${PROTEINS_FAA}. Run init first."; exit 1; }

    diamond blastp \
      --query "${PROTEINS_FAA}" \
      --db "${SWISSPROT_DB}" \
      --out "${ORTHO_TSV}" \
      --outfmt 6 qseqid sseqid stitle pident length qlen slen qstart qend sstart send evalue bitscore \
      --max-target-seqs 1 \
      --evalue 1e-5 \
      --threads "${THREADS}" \
      --sensitive

    python workflow/update_global_annotation_table.py merge-orthology \
      --table-tsv "${TABLE_TSV}" \
      --diamond-tsv "${ORTHO_TSV}"
    ;;

  domains)
    [[ -s "${PROTEINS_FAA}" ]] || { echo "ERROR: proteins FASTA missing: ${PROTEINS_FAA}. Run init first."; exit 1; }

    hmmscan \
      --cpu "${THREADS}" \
      --domtblout "${DOMAIN_TBLOUT}" \
      "${PFAM_DB}" \
      "${PROTEINS_FAA}" \
      > "${OUTDIR}/03_domains/hmmscan.stdout.txt"

    python workflow/update_global_annotation_table.py merge-domains \
      --table-tsv "${TABLE_TSV}" \
      --domtblout "${DOMAIN_TBLOUT}"
    ;;

  functional)
    [[ -s "${PROTEINS_FAA}" ]] || { echo "ERROR: proteins FASTA missing: ${PROTEINS_FAA}. Run init first."; exit 1; }

    emapper.py \
      -i "${PROTEINS_FAA}" \
      --itype proteins \
      --cpu "${THREADS}" \
      --data_dir "${EGGNOG_DATA_DIR}" \
      --output "$(basename "${EGGNOG_PREFIX}")" \
      --output_dir "$(dirname "${EGGNOG_PREFIX}")" \
      --override

    python workflow/update_global_annotation_table.py merge-functional \
      --table-tsv "${TABLE_TSV}" \
      --eggnog-annotations "${EGGNOG_PREFIX}.emapper.annotations"
    ;;

  taxonomy)
    python workflow/update_global_annotation_table.py merge-taxonomy \
      --table-tsv "${TABLE_TSV}" \
      --diamond-tsv "${ORTHO_TSV}" \
      --eggnog-annotations "${EGGNOG_PREFIX}.emapper.annotations"
    ;;

  *)
    echo "ERROR: unsupported step '${STEP}'"
    exit 1
    ;;
esac

echo "Finished step: ${STEP}"
