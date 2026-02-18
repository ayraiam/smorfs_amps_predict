#!/usr/bin/env bash
# ==========================================================
# Script: workflow/downstream_analysis.sh
# Purpose:
#   - Create (once) a conda env for downstream processing
#   - Step A (Python): parse MetaFlye flye.log -> metrics TSV
#   - Step B (R): boxplots for each metrics column
#   - Step C (Bash): build global nonredundant peptide set + AMP prediction hook
#   - Step D (Bash): build global nonredundant contig set + map reads back per replicate
#
# Default inputs:
#   results/assembly_metaflye/*/flye.log
#   results/smorfs/*/catalog/nr_catalog_rep_seq.fasta
#   results/assembly_metaflye/*/assembly.fasta
#   metadata/metagenome_files.txt
#
# Default outputs:
#   results/assembly_metaflye/finalize_metrics.tsv
#   results/assembly_metaflye/finalize_boxplots/boxplot_*.png
#   results/catalog_global/global_nr_peptides.faa (+ clusters)
#   results/catalog_global/global_nr_contigs.fna (+ mapping tables)
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

# New global catalog outputs
GLOBAL_DIR="${GLOBAL_DIR:-${BASE_RESULTS_DIR}/catalog_global}"
GLOBAL_PEPTIDES="${GLOBAL_PEPTIDES:-${GLOBAL_DIR}/global_nr_peptides.faa}"
GLOBAL_CONTIGS="${GLOBAL_CONTIGS:-${GLOBAL_DIR}/global_nr_contigs.fna}"
MAP_OUTDIR="${MAP_OUTDIR:-${GLOBAL_DIR}/mapping}"

# Metadata map for replicates
METADATA_MAP="${METADATA_MAP:-metadata/metagenome_files.txt}"
SAMPLE_COL="${SAMPLE_COL:-SampleID}"
FASTQ_COL="${FASTQ_COL:-FASTQ_Filename}"

# Run flags
RUN_PARSE=1
RUN_PLOTS=1
RUN_AMPS=0
RUN_MAP=0
RUN_MACREL_ATTACH=0

# Macrel attach (per-sample predicted_smorfs.tsv)
RUN_MACREL_ATTACH=0
PREDICTED_SMORFS_TSV=""
MACREL_ID_COL="feature_id"
MACREL_SEQ_COL=""
MACREL_ATTACH_OUT=""

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

    # Python for parsing + R for plotting + tools for global NR + mapping
    "${installer}" create -y -p "${ENV_PREFIX}" --override-channels \
      -c conda-forge -c bioconda \
      python=3.10 \
      pandas \
      matplotlib \
      r-base=4.3 \
      r-ggplot2 \
      r-readr \
      r-dplyr \
      r-ggbeeswarm \
      mmseqs2 \
      minimap2 \
      samtools \
      seqkit \
      bedtools \
      macrel

    msg "Env created: ${ENV_PREFIX}"

    trap - EXIT INT TERM
    cleanup_lock
  fi

  set +u
  conda activate "${ENV_PREFIX}"
  set -u

  msg "Python: $(python -V)"
  msg "R:      $(R --version | head -n 1)"
  msg "mmseqs: $(mmseqs version 2>/dev/null || echo '<ok>')"
  msg "minimap2: $(minimap2 --version 2>/dev/null || echo '<ok>')"
  msg "samtools: $(samtools --version | head -n1)"
}

# ---------------------------
# Step A/B
# ---------------------------
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

# ---------------------------
# Step C: Global NR peptides + AMP prediction hook
# ---------------------------
step_c_amps() {
  msg "Step C: build global NR peptide set + AMP prediction hook"
  mkdir -p "${GLOBAL_DIR}"

  # Gather per-site peptide NR catalogs
  local tmp_all="${GLOBAL_DIR}/_all_site_nr_peptides.faa"
  : > "${tmp_all}"

  shopt -s nullglob
  local files=( "${BASE_RESULTS_DIR}"/smorfs/*/catalog/nr_catalog_rep_seq.fasta )
  shopt -u nullglob

  [[ ${#files[@]} -gt 0 ]] || die "No peptide catalogs found at ${BASE_RESULTS_DIR}/smorfs/*/catalog/nr_catalog_rep_seq.fasta"

  msg "Found ${#files[@]} per-site peptide catalogs. Concatenating..."
  cat "${files[@]}" > "${tmp_all}"

  # Global NR with mmseqs (protein clustering)
  local mmdb="${GLOBAL_DIR}/mmseqs_pep_db"
  local tmpdir="${GLOBAL_DIR}/mmseqs_tmp"
  local clust="${GLOBAL_DIR}/mmseqs_pep_clusters"
  mkdir -p "${tmpdir}"

  msg "Running mmseqs easy-cluster (proteins) -> global NR representatives"
  # Conservative defaults; you can tweak --min-seq-id and -c later
  mmseqs easy-cluster "${tmp_all}" "${clust}" "${tmpdir}" --min-seq-id 0.95 -c 0.8 >/dev/null

  # mmseqs easy-cluster writes reps in:
  #   ${clust}_rep_seq.fasta
  [[ -f "${clust}_rep_seq.fasta" ]] || die "Expected reps fasta not found: ${clust}_rep_seq.fasta"

  cp "${clust}_rep_seq.fasta" "${GLOBAL_PEPTIDES}"
  msg "Wrote global NR peptides: ${GLOBAL_PEPTIDES}"

  # AMP prediction hook (optional): if you add your own predictor script, we call it.
  # You can implement workflow/predict_amps.py later; for now we just log.
  msg "Running Macrel AMP prediction on global NR peptides"
  python workflow/predict_amps.py "${GLOBAL_PEPTIDES}" \
    --out "${GLOBAL_DIR}/amp_predictions.tsv" \
    --threads "${CPUS:-8}"
  msg "Wrote AMP predictions: ${GLOBAL_DIR}/amp_predictions.tsv"

}

step_c_macrel_attach_predicted_smorfs() {
  [[ -n "${PREDICTED_SMORFS_TSV}" ]] || die "--predicted-smorfs not provided"
  [[ -f "${PREDICTED_SMORFS_TSV}" ]] || die "Missing: ${PREDICTED_SMORFS_TSV}"

  local outdir
  outdir="$(dirname "${PREDICTED_SMORFS_TSV}")"

  local peptides_faa="${outdir}/predicted_smorfs.peptides_for_macrel.faa"
  local rejected_tsv="${outdir}/predicted_smorfs.peptides_for_macrel.rejected.tsv"
  local macrel_tsv="${outdir}/macrel_predictions.normalized.tsv predominan"
  local merged_tsv

  if [[ -n "${MACREL_ATTACH_OUT:-}" ]]; then
    merged_tsv="${MACREL_ATTACH_OUT}"
  else
    merged_tsv="${outdir}/predicted_smorfs.with_macrel.tsv"
  fi

  msg "Building peptide FASTA for Macrel from: ${PREDICTED_SMORFS_TSV}"
  msg "  ID col : ${MACREL_ID_COL:-feature_id}"
  msg "  SEQ col: ${MACREL_SEQ_COL}"
  msg "  FASTA  : ${peptides_faa}"
  msg "  Reject : ${rejected_tsv}"

  python - "${PREDICTED_SMORFS_TSV}" "${peptides_faa}" "${rejected_tsv}" "${MACREL_ID_COL:-feature_id}" "${MACREL_SEQ_COL}" <<'PY'
import sys
import pandas as pd
from pathlib import Path

tsv = Path(sys.argv[1])
faa = Path(sys.argv[2])
rej = Path(sys.argv[3])
id_col = sys.argv[4]
seq_col = sys.argv[5]

df = pd.read_csv(tsv, sep="\t", dtype=str)

if id_col not in df.columns:
    raise SystemExit(f"ERROR: id-col '{id_col}' not found. Columns={df.columns.tolist()}")
if seq_col not in df.columns:
    raise SystemExit(f"ERROR: seq-col '{seq_col}' not found. Columns={df.columns.tolist()}")

df[id_col] = df[id_col].astype(str).str.strip()
df[seq_col] = df[seq_col].astype(str).str.strip()

if df[id_col].duplicated().any():
    dups = df.loc[df[id_col].duplicated(), id_col].head(10).tolist()
    raise SystemExit(f"ERROR: duplicated IDs in {id_col}. Examples: {dups}")

ok = set("ACDEFGHIKLMNPQRSTVWY")

n_out = 0
n_bad = 0
n_empty = 0

with faa.open("w") as out, rej.open("w") as r:
    r.write("peptide_id\tbad_chars\tpeptide_seq\n")

    for pid, seq in zip(df[id_col], df[seq_col]):
        if not pid or str(pid).lower() == "nan":
            n_empty += 1
            continue
        if not seq or str(seq).lower() == "nan":
            n_empty += 1
            continue

        pid = str(pid).strip()
        seq = str(seq).strip().replace("\r", "").replace(" ", "").upper().rstrip("*")

        if not seq:
            n_empty += 1
            continue

        bad = sorted({c for c in seq if c not in ok})
        if bad:
            n_bad += 1
            r.write(f"{pid}\t{''.join(bad)}\t{seq}\n")
            continue

        out.write(f">{pid}\n{seq}\n")
        n_out += 1

print(f"Wrote FASTA: {faa} (kept={n_out}, rejected={n_bad}, empty={n_empty})", file=sys.stderr)
print(f"Rejected list: {rej}", file=sys.stderr)
PY

  msg "Running Macrel via workflow/predict_amps.py"
  python workflow/predict_amps.py "${peptides_faa}" \
    --out "${macrel_tsv}" \
    --threads "${CPUS:-8}"

  msg "Attaching Macrel predictions back to predicted_smorfs.tsv"
  python workflow/attach_macrel_to_predicted_smorfs.py \
    --predicted "${PREDICTED_SMORFS_TSV}" \
    --macrel "${macrel_tsv}" \
    --out "${merged_tsv}" \
    --id-col "${MACREL_ID_COL:-feature_id}" \
    --seq-col "${MACREL_SEQ_COL}"

  msg "Wrote merged file: ${merged_tsv}"
  msg "If Macrel still fails, inspect rejected peptides: ${rejected_tsv}"
}

# ---------------------------
# Step D: Global NR contigs + mapping per replicate
# ---------------------------
step_d_map() {
  msg "Step D: build global NR contigs + map reads back per replicate"
  mkdir -p "${GLOBAL_DIR}" "${MAP_OUTDIR}"

  # Gather assemblies
  shopt -s nullglob
  local asms=( "${BASE_RESULTS_DIR}"/assembly_metaflye/*/assembly.fasta )
  shopt -u nullglob
  [[ ${#asms[@]} -gt 0 ]] || die "No assemblies found at ${BASE_RESULTS_DIR}/assembly_metaflye/*/assembly.fasta"

  local all_contigs="${GLOBAL_DIR}/_all_site_contigs.fna"
  cat "${asms[@]}" > "${all_contigs}"

  # Deduplicate contigs with mmseqs nucleotide clustering (search-type 3)
  local tmpdir="${GLOBAL_DIR}/mmseqs_tmp_nt"
  local clust="${GLOBAL_DIR}/mmseqs_nt_clusters"
  mkdir -p "${tmpdir}"

  msg "Running mmseqs easy-cluster (nucleotide contigs) -> global NR contigs"
  mmseqs easy-cluster "${all_contigs}" "${clust}" "${tmpdir}" --search-type 3 --min-seq-id 0.99 -c 0.9 >/dev/null
  [[ -f "${clust}_rep_seq.fasta" ]] || die "Expected reps fasta not found: ${clust}_rep_seq.fasta"

  cp "${clust}_rep_seq.fasta" "${GLOBAL_CONTIGS}"
  msg "Wrote global NR contigs: ${GLOBAL_CONTIGS}"

  # Build minimap2 index (optional; minimap2 can map without explicit index, but this helps)
  local idx="${GLOBAL_CONTIGS}.mmi"
  if [[ ! -f "${idx}" ]]; then
    msg "Indexing global NR contigs with minimap2"
    minimap2 -d "${idx}" "${GLOBAL_CONTIGS}" >/dev/null
  fi

  # Metadata map required for per-replicate fastqs
  [[ -f "${METADATA_MAP}" ]] || die "Missing metadata map: ${METADATA_MAP}"

  # Determine FASTQ_DIR if provided; default data/
  local fastq_dir="${FASTQ_DIR:-data}"
  [[ -d "${fastq_dir}" ]] || die "FASTQ_DIR not found: ${fastq_dir}"

  # Output summary tables per replicate:
  #   idxstats: contig lengths + mapped reads
  #   depth: mean depth per contig
  msg "Mapping each replicate separately (per row in ${METADATA_MAP})"
  tail -n +2 "${METADATA_MAP}" | while IFS=$'\t' read -r line; do
    # Robust column picking: assume tab-delimited with headers; simplest parse via awk for column indices
    break
  done

  # Column indices
  local sidx fidx
  sidx=$(awk -v c="${SAMPLE_COL}" 'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i==c) print i}' "${METADATA_MAP}")
  fidx=$(awk -v c="${FASTQ_COL}"  'BEGIN{FS="\t"} NR==1{for(i=1;i<=NF;i++) if($i==c) print i}' "${METADATA_MAP}")
  [[ -n "${sidx}" ]] || die "Could not find SAMPLE_COL=${SAMPLE_COL} in header of ${METADATA_MAP}"
  [[ -n "${fidx}" ]] || die "Could not find FASTQ_COL=${FASTQ_COL} in header of ${METADATA_MAP}"

  tail -n +2 "${METADATA_MAP}" | awk -v FS="\t" -v s="${sidx}" -v f="${fidx}" '{print $s "\t" $f}' | while IFS=$'\t' read -r sample fastq; do
    [[ -n "${sample}" ]] || continue
    [[ -n "${fastq}" ]] || continue

    local fq="${fastq_dir}/${fastq}"
    [[ -f "${fq}" ]] || { msg "WARN: missing FASTQ ${fq}; skipping"; continue; }

    local outprefix="${MAP_OUTDIR}/${sample}"
    local bam="${outprefix}.bam"
    local sorted="${outprefix}.sorted.bam"
    local idxstats="${outprefix}.idxstats.tsv"
    local depth="${outprefix}.depth.tsv"
    local meandepth="${outprefix}.mean_depth_per_contig.tsv"

    if [[ -f "${sorted}" ]]; then
      msg "[${sample}] sorted BAM exists; skipping mapping"
    else
      msg "[${sample}] minimap2 mapping -> BAM"
      minimap2 -t "${CPUS:-8}" -ax map-ont "${idx}" "${fq}" \
        | samtools view -@ "${CPUS:-8}" -b -o "${bam}" -

      samtools sort -@ "${CPUS:-8}" -o "${sorted}" "${bam}"
      samtools index "${sorted}"
      rm -f "${bam}"
    fi

    msg "[${sample}] idxstats"
    samtools idxstats "${sorted}" > "${idxstats}"

    msg "[${sample}] depth"
    samtools depth -aa "${sorted}" > "${depth}"

    # Mean depth per contig
    awk 'BEGIN{OFS="\t"} {sum[$1]+=$3; n[$1]+=1} END{for (c in sum) print c, sum[c]/n[c]}' "${depth}" \
      | sort -k1,1 > "${meandepth}"

    msg "[${sample}] wrote:"
    msg "  ${idxstats}"
    msg "  ${meandepth}"
  done

  msg "Mapping outputs in: ${MAP_OUTDIR}"
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

  --parse-only           Run only Step A
  --plots-only           Run only Step B
  --amps-only            Run only Step C
  --map-only             Run only Step D

  --run-amps             Enable Step C in addition to A/B
  --run-map              Enable Step D in addition to A/B

  --global-dir PATH      Global output dir (default: results/catalog_global)
  --metadata-map PATH    Replicate mapping file (default: metadata/metagenome_files.txt)
  --sample-col STR       SampleID column name in metadata map (default: SampleID)
  --fastq-col STR        FASTQ filename column name (default: FASTQ_Filename)

  --macrel-attach-only       Run ONLY: Macrel on predicted_smorfs.tsv and attach results
  --predicted-smorfs PATH    Input predicted smORFs TSV
  --id-col STR               ID column in predicted TSV (recommended: feature_id)
  --seq-col STR              AA sequence column in predicted TSV (e.g. aa_seq)
  --macrel-attach-out PATH   Output merged TSV (default: predicted_smorfs.with_macrel.tsv)

  -h, --help             Show help

Notes:
  - Step C expects: results/smorfs/*/catalog/nr_catalog_rep_seq.fasta
  - Step D expects: results/assembly_metaflye/*/assembly.fasta and metadata/metagenome_files.txt
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

    --global-dir) GLOBAL_DIR="$2"; GLOBAL_PEPTIDES="${GLOBAL_DIR}/global_nr_peptides.faa"; GLOBAL_CONTIGS="${GLOBAL_DIR}/global_nr_contigs.fna"; MAP_OUTDIR="${GLOBAL_DIR}/mapping"; shift 2 ;;
    --metadata-map) METADATA_MAP="$2"; shift 2 ;;
    --sample-col) SAMPLE_COL="$2"; shift 2 ;;
    --fastq-col) FASTQ_COL="$2"; shift 2 ;;

    --parse-only) RUN_PLOTS=0; RUN_AMPS=0; RUN_MAP=0; shift ;;
    --plots-only) RUN_PARSE=0; RUN_AMPS=0; RUN_MAP=0; shift ;;
    --amps-only) RUN_PARSE=0; RUN_PLOTS=0; RUN_MAP=0; RUN_AMPS=1; shift ;;
    --map-only) RUN_PARSE=0; RUN_PLOTS=0; RUN_AMPS=0; RUN_MAP=1; shift ;;

    --run-amps) RUN_AMPS=1; shift ;;
    --run-map) RUN_MAP=1; shift ;;

    --predicted-smorfs) PREDICTED_SMORFS_TSV="$2"; shift 2 ;;
    --macrel-out) MACREL_ATTACH_OUT="$2"; shift 2 ;;

    --macrel-attach-only)
    RUN_PARSE=0
    RUN_PLOTS=0
    RUN_AMPS=0
    RUN_MAP=0
    RUN_MACREL_ATTACH=1
    shift
    ;;
    --predicted-smorfs) PREDICTED_SMORFS_TSV="$2"; shift 2 ;;
    --id-col) MACREL_ID_COL="$2"; shift 2 ;;
    --seq-col) MACREL_SEQ_COL="$2"; shift 2 ;;
    --macrel-attach-out) MACREL_ATTACH_OUT="$2"; shift 2 ;;

    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1 (use --help)" ;;
  esac
done

if [[ "${RUN_MACREL_ATTACH}" -eq 1 ]]; then
  [[ -n "${PREDICTED_SMORFS_TSV}" ]] || die "--predicted-smorfs is required with --macrel-attach-only"
  [[ -n "${MACREL_SEQ_COL}" ]] || die "--seq-col is required with --macrel-attach-only"
fi

# ---------------------------
# Main
# ---------------------------
ensure_env_once

[[ -f workflow/summarize_flye_logs.py ]] || die "Missing: workflow/summarize_flye_logs.py"
[[ -f workflow/plot_metaflye_metrics.R ]] || die "Missing: workflow/plot_metaflye_metrics.R"
[[ -f workflow/predict_amps.py ]] || die "Missing: workflow/predict_amps.py"
[[ -f workflow/attach_macrel_to_predicted_smorfs.py ]] || die "Missing: workflow/attach_macrel_to_predicted_smorfs.py"

if [[ "${RUN_PARSE}" -eq 1 ]]; then step_a_parse; fi
if [[ "${RUN_PLOTS}" -eq 1 ]]; then step_b_plots_r; fi
if [[ "${RUN_AMPS}" -eq 1 ]]; then step_c_amps; fi
if [[ "${RUN_MACREL_ATTACH}" -eq 1 ]]; then step_c_macrel_attach_predicted_smorfs; fi
if [[ "${RUN_MAP}" -eq 1 ]]; then step_d_map; fi

msg "DONE"
