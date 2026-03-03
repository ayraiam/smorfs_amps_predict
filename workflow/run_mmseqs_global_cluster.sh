#!/usr/bin/env bash
set -euo pipefail

# ==========================================================
# Script: workflow/run_mmseqs_global_cluster.sh
# Purpose:
#   1) Prefix peptide FASTA headers with environment label to avoid ID collisions
#   2) Concatenate into one global FASTA
#   3) Run MMseqs2 clustering
#   4) Export cluster assignments as TSVs
#
# Outputs (default):
#   results/smorfs/GLOBAL/catalog/peptides_all_global.faa
#   results/smorfs/GLOBAL/mmseqs/clusters_pairs.tsv   (rep <tab> member)
#   results/smorfs/GLOBAL/mmseqs/cluster_map.tsv      (feature_id <tab> cluster_id)
# ==========================================================

# Defaults
RESULTS_DIR="${RESULTS_DIR:-results}"
OUT_BASE="${RESULTS_DIR}/smorfs/GLOBAL"
CATALOG_DIR="${OUT_BASE}/catalog"
MMSEQS_DIR="${OUT_BASE}/mmseqs"
TMP_DIR="${MMSEQS_DIR}/tmp"

# Environments (edit here if you rename)
ENV_NAMES=("RIPARIA" "PENEIRA" "CAMPINA" "UNIFORME")

# MMseqs params (sane defaults; adjust later)
MIN_SEQ_ID="${MIN_SEQ_ID:-0.95}"
COV="${COV:-0.8}"
COV_MODE="${COV_MODE:-1}"
THREADS="${THREADS:-8}"

die(){ echo "ERROR: $*" >&2; exit 1; }
msg(){ echo "[$(date +'%F %T')] $*" >&2; }

usage(){
  cat <<EOF
Usage:
  bash workflow/run_mmseqs_global_cluster.sh [options]

Options:
  --results-dir DIR     (default: results)
  --min-seq-id FLOAT    (default: 0.95)
  --cov FLOAT           (default: 0.8)
  --cov-mode INT        (default: 1)
  --threads INT         (default: 8)

Notes:
  Expects input FASTAs at:
    results/smorfs/<ENV>_GLOBAL/catalog/peptides_all.faa
  and writes a global FASTA + mmseqs outputs under:
    results/smorfs/GLOBAL/
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --min-seq-id) MIN_SEQ_ID="$2"; shift 2 ;;
    --cov) COV="$2"; shift 2 ;;
    --cov-mode) COV_MODE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown arg: $1" ;;
  esac
done

# Paths (recompute if RESULTS_DIR overridden)
OUT_BASE="${RESULTS_DIR}/smorfs/GLOBAL"
CATALOG_DIR="${OUT_BASE}/catalog"
MMSEQS_DIR="${OUT_BASE}/mmseqs"
TMP_DIR="${MMSEQS_DIR}/tmp"

mkdir -p "${CATALOG_DIR}" "${MMSEQS_DIR}" "${TMP_DIR}"

GLOBAL_FASTA="${CATALOG_DIR}/peptides_all_global.faa"

msg "============================================"
msg "Global clustering: building prefixed FASTAs"
msg "RESULTS_DIR   : ${RESULTS_DIR}"
msg "OUT_BASE      : ${OUT_BASE}"
msg "GLOBAL_FASTA  : ${GLOBAL_FASTA}"
msg "MMseqs params : --min-seq-id ${MIN_SEQ_ID} -c ${COV} --cov-mode ${COV_MODE}"
msg "THREADS       : ${THREADS}"
msg "============================================"

# 1) Prefix and concatenate
: > "${GLOBAL_FASTA}"  # truncate

for ENV in "${ENV_NAMES[@]}"; do
  IN_FA="${RESULTS_DIR}/smorfs/${ENV}_GLOBAL/catalog/peptides_all.faa"
  [[ -f "${IN_FA}" ]] || die "Missing input FASTA: ${IN_FA}"

  msg "Prefixing: ${ENV}  <- ${IN_FA}"
  # Prefix only the first token after '>' (keeps metadata after spaces intact)
  awk -v env="${ENV}" '
    /^>/ {
      sub(/^>/, ">");
      # split header into first token and the rest (preserve rest)
      split($0, a, " ");
      id = substr(a[1], 2);  # drop leading >
      rest = "";
      # rebuild rest of line exactly
      if (length($0) > length(a[1])) {
        rest = substr($0, length(a[1]) + 1);
      }
      print ">" env "|" id rest;
      next
    }
    { print }
  ' "${IN_FA}" >> "${GLOBAL_FASTA}"
done

# 2) Verify no duplicate FASTA IDs
msg "Checking for duplicate FASTA IDs after prefixing..."
DUPES=$(grep '^>' "${GLOBAL_FASTA}" | cut -d' ' -f1 | sort | uniq -d | head -n 20 || true)
if [[ -n "${DUPES}" ]]; then
  echo "${DUPES}" >&2
  die "Duplicate IDs still present after prefixing. Investigate above."
fi
msg "OK: no duplicate FASTA IDs detected."

# 3) MMseqs2 clustering
command -v mmseqs >/dev/null 2>&1 || die "mmseqs not found in PATH. Activate your env that has mmseqs2."

DB="${MMSEQS_DIR}/peptides_db"
CLU="${MMSEQS_DIR}/clusters"
PAIRS_TSV="${MMSEQS_DIR}/clusters_pairs.tsv"   # rep <tab> member
MAP_TSV="${MMSEQS_DIR}/cluster_map.tsv"        # feature_id <tab> cluster_id

msg "Running mmseqs createdb..."
mmseqs createdb "${GLOBAL_FASTA}" "${DB}"

msg "Running mmseqs cluster..."
mmseqs cluster "${DB}" "${CLU}" "${TMP_DIR}" \
  --min-seq-id "${MIN_SEQ_ID}" \
  -c "${COV}" \
  --cov-mode "${COV_MODE}" \
  --threads "${THREADS}"

msg "Exporting rep-member pairs TSV..."
mmseqs createtsv "${DB}" "${DB}" "${CLU}" "${PAIRS_TSV}"

msg "Building feature_id -> cluster_id map TSV..."
# PAIRS_TSV: rep \t member
# We want:   feature_id(member) \t cluster_id(rep)
awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1}' "${PAIRS_TSV}" > "${MAP_TSV}"

msg "DONE."
msg "Outputs:"
msg "  ${GLOBAL_FASTA}"
msg "  ${PAIRS_TSV}"
msg "  ${MAP_TSV}"
