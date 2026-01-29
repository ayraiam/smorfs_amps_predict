#!/usr/bin/env bash
# ==========================================================
# Script: workflow/run_smorfs_pipeline.sh
# Purpose:
#   - Create env with tools for bacterial + fungal smORF discovery
#   - For each MetaFlye assembly.fasta:
#       A) Flag contigs <500 bp as low-confidence (DO NOT remove)
#       B) Track per-site assembly stats in single TSV (append-safe)
#       C) Classify contigs (Tiara) -> branch bac vs fungal
#       D) Prodigal meta-mode ORFs for bac contigs
#       E) SmORFinder for bac contigs
#       F) Funannotate for fungal contigs + flag <=100 aa peptides
#       G) Pool peptides and de-replicate (MMseqs2 easy-cluster)
#
# Inputs:
#   results/assembly_metaflye/<sample_name>/assembly.fasta
#
# Outputs (per sample):
#   results/smorfs/<sample_name>/
#     contigs/ (split + classified)
#     bac/ (prodigal + smorfinder)
#     fungi/ (funannotate)
#     catalog/ (pooled + nonredundant clusters)
#
# Global outputs:
#   results/smorfs/assembly_stats.tsv         (one row per sample)
#   results/smorfs/contig_stats.tsv           (one row per contig)
# ==========================================================
set -euo pipefail

# ---------------------------
# Defaults
# ---------------------------
# Env naming/prefix (HPC-safe, like MetaFlye)
ENV_NAME="${SMORFS_ENV:-smorfs_amps_env}"
ENV_PREFIX_DIR="${ENV_PREFIX_DIR:-envs}"
ENV_PREFIX="${ENV_PREFIX_DIR}/${ENV_NAME}"

# IMPORTANT:
# RESULTS_DIR passed from runall.sh is the BASE results folder (default: "results")
BASE_RESULTS_DIR="${RESULTS_DIR:-results}"

# Where smORFs outputs go:
RESULTS_DIR="${BASE_RESULTS_DIR}/smorfs"

# Where MetaFlye assemblies are:
ASSEMBLY_ROOT="${BASE_RESULTS_DIR}/assembly_metaflye"

MIN_CONTIG_BP=500
MAX_FUNGAL_PEPTIDE_AA=100
CPUS=8

# For funannotate DB (optional but strongly recommended)
# If empty, funannotate may still run but often complains.
FUNANNOTATE_DB_DIR="${FUNANNOTATE_DB_DIR:-}"

# ---------------------------
# Helpers
# ---------------------------
die() { echo "ERROR: $*" >&2; exit 1; }
msg() { echo "[`date +'%F %T'`] $*" >&2; }

need_cmd() {
  command -v "$1" >/dev/null 2>&1 || die "Missing required command in PATH: $1"
}

mkdirp() { mkdir -p "$1"; }

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
    # Acquire lock
    while ! mkdir "${lockdir}" 2>/dev/null; do
      msg "Waiting for env lock..."
      sleep 5
    done

    # Another process may have created it while we waited
    if [[ -d "${ENV_PREFIX}" ]]; then
      msg "Env appeared while waiting. Continuing."
      rmdir "${lockdir}" || true
    else
      msg "Creating conda env at: ${ENV_PREFIX}"
      local installer="conda"
      if have_cmd mamba; then installer="mamba"; fi

      # Use strict channel priority to reduce solver pain
      "${installer}" config set channel_priority strict >/dev/null 2>&1 || true

      # Create a PREFIX env (NOT named env) like MetaFlye
      "${installer}" create -y -p "${ENV_PREFIX}" --override-channels \
        -c conda-forge -c bioconda \
        python=3.8 \
        prodigal \
        mmseqs2 \
        seqkit \
        csvtk \
        pigz \
        tiara

      msg "Env created: ${ENV_PREFIX}"
      rmdir "${lockdir}" || true
    fi
  fi

  # conda activation scripts may reference unset vars (e.g., MKL_INTERFACE_LAYER)
  # so temporarily disable nounset just for activation
  set +u
  conda activate "${ENV_PREFIX}"
  set -u

  # Sanity checks
  have_cmd seqkit   || die "seqkit not found after activating env."
  have_cmd prodigal || die "prodigal not found after activating env."
  have_cmd mmseqs   || die "mmseqs not found after activating env."
  have_cmd tiara    || die "tiara not found after activating env."
  have_cmd pigz     || die "pigz not found after activating env."
}

# ---------------------------
# 1) Create env
# ---------------------------
create_env() {
  msg "Ensuring smORFs env exists (prefix mode): ${ENV_PREFIX}"
  ensure_env_once

  msg "Installing pip tools (SmORFinder + funannotate) into env..."
  python -m pip install --upgrade pip
  python -m pip install smorfinder funannotate

  # Confirm CLIs are available (best-effort)
  command -v smorf >/dev/null 2>&1 || msg "WARNING: 'smorf' not found after pip install. Check pip output."
  command -v funannotate >/dev/null 2>&1 || msg "WARNING: 'funannotate' not found after pip install. Check pip output."

  msg "Env ready. To activate:"
  msg "  source \"\$(conda info --base)/etc/profile.d/conda.sh\""
  msg "  conda activate \"${ENV_PREFIX}\""
  msg "Optional (recommended) funannotate DB setup:"
  msg "  funannotate setup -d /path/to/funannotate_db   (set FUNANNOTATE_DB_DIR=/path/to/funannotate_db)"
}

# ---------------------------
# 2) Flag low-confidence contigs + track stats TSVs
# ---------------------------
init_tables() {
  mkdirp "${RESULTS_DIR}"
  if [[ ! -f "${RESULTS_DIR}/assembly_stats.tsv" ]]; then
    printf "sample_id\tassembly_fasta\tcontigs_total\ttotal_bp\tmin_len\tmax_len\tn50\tcontigs_lt_%dbp\n" \
      "${MIN_CONTIG_BP}" > "${RESULTS_DIR}/assembly_stats.tsv"
  fi
  if [[ ! -f "${RESULTS_DIR}/contig_stats.tsv" ]]; then
    printf "sample_id\tcontig_id\tlength_bp\tlow_confidence\n" > "${RESULTS_DIR}/contig_stats.tsv"
  fi
}

flag_contigs_and_record_stats() {
  local sample_id="$1"
  local asm_fa="$2"
  local outdir="$3"

  mkdirp "${outdir}/contigs"

  # Per-contig lengths (seqkit) -> contig_stats.tsv
  # Also create list of low-confidence contigs (<MIN_CONTIG_BP)
  msg "[${sample_id}] Computing contig lengths + low-confidence flags (threshold=${MIN_CONTIG_BP} bp)"
  seqkit fx2tab -nl "${asm_fa}" \
    | awk -v S="${sample_id}" -v MIN="${MIN_CONTIG_BP}" 'BEGIN{OFS="\t"}{
        id=$1; len=$2; low=(len<MIN?"TRUE":"FALSE");
        print S,id,len,low
      }' >> "${RESULTS_DIR}/contig_stats.tsv"

  seqkit fx2tab -nl "${asm_fa}" \
    | awk -v MIN="${MIN_CONTIG_BP}" '$2<MIN{print $1}' > "${outdir}/contigs/low_confidence.ids"

  # Assembly summary stats (seqkit stats -T gives a TSV line)
  # We'll parse out key fields.
  msg "[${sample_id}] Computing assembly summary stats"
  # seqkit stats -T output columns include: file, format, type, num_seqs, sum_len, min_len, avg_len, max_len, Q1,Q2,Q3,sum_gap,N50,Q20,Q30
  local stats_tsv
  stats_tsv="$(seqkit stats -T "${asm_fa}" | tail -n 1)"

  local contigs_total total_bp min_len max_len n50
  contigs_total="$(echo "${stats_tsv}" | awk -F'\t' '{print $4}')"
  total_bp="$(echo "${stats_tsv}" | awk -F'\t' '{print $5}')"
  min_len="$(echo "${stats_tsv}" | awk -F'\t' '{print $6}')"
  max_len="$(echo "${stats_tsv}" | awk -F'\t' '{print $9}')"
  n50="$(echo "${stats_tsv}" | awk -F'\t' '{print $14}')"

  local lt_min
  lt_min="$(seqkit fx2tab -nl "${asm_fa}" | awk -v MIN="${MIN_CONTIG_BP}" '$2<MIN{c++} END{print c+0}')"

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${sample_id}" "${asm_fa}" "${contigs_total}" "${total_bp}" "${min_len}" "${max_len}" "${n50}" "${lt_min}" \
    >> "${RESULTS_DIR}/assembly_stats.tsv"
}

# ---------------------------
# 3) Classify contigs (flag likely fungal, then branch)
# ---------------------------
# Strategy:
#   Use Tiara to classify each contig: Bacteria/Archaea/Prokarya/Eukarya/Organelle/Unknown. :contentReference[oaicite:9]{index=9}
#   Then:
#     - "bacterial bucket": Bacteria + Archaea + Prokarya
#     - "fungal bucket":   Eukarya
#     - "other bucket":    Organelle + Unknown (kept, not discarded)
#
# Why Tiara?
#   It's designed for assembled contigs and separates eukaryotic nuclear sequences from prokaryotes,
#   with special handling of organelles. :contentReference[oaicite:10]{index=10}
classify_contigs_tiara() {
  local sample_id="$1"
  local asm_fa="$2"
  local outdir="$3"

  mkdirp "${outdir}/contigs/classify"
  msg "[${sample_id}] Classifying contigs with Tiara"
  # Tiara CLI may differ slightly by version; commonly:
  #   tiara -i assembly.fasta -o tiara.tsv -t <threads>
  # We'll use a conservative call and fail loudly if the CLI differs.
  tiara -i "${asm_fa}" -o "${outdir}/contigs/classify/tiara.tsv" -t "${CPUS}"

  # tiara.tsv typically: contig_id <tab> label
  # Split contigs by label -> FASTA files
  msg "[${sample_id}] Splitting contigs by Tiara label"
  awk -F'\t' 'BEGIN{OFS="\t"} NR>0{print $1,$2}' "${outdir}/contigs/classify/tiara.tsv" \
    > "${outdir}/contigs/classify/tiara.map.tsv"

  awk -F'\t' '$2 ~ /Eukarya/i {print $1}' "${outdir}/contigs/classify/tiara.map.tsv" > "${outdir}/contigs/fungi.ids"
  awk -F'\t' '$2 ~ /(Bacteria|Archaea|Prokarya)/i {print $1}' "${outdir}/contigs/classify/tiara.map.tsv" > "${outdir}/contigs/bac.ids"
  awk -F'\t' '$2 ~ /(Organelle|Unknown)/i {print $1}' "${outdir}/contigs/classify/tiara.map.tsv" > "${outdir}/contigs/other.ids"

  # Extract FASTAs
  seqkit grep -f "${outdir}/contigs/bac.ids"   "${asm_fa}" > "${outdir}/contigs/bac_contigs.fasta"   || true
  seqkit grep -f "${outdir}/contigs/fungi.ids" "${asm_fa}" > "${outdir}/contigs/fungi_contigs.fasta" || true
  seqkit grep -f "${outdir}/contigs/other.ids" "${asm_fa}" > "${outdir}/contigs/other_contigs.fasta" || true

  # Record basic counts in a small note file
  {
    echo "sample_id=${sample_id}"
    echo "bac_contigs=$(grep -c '^>' "${outdir}/contigs/bac_contigs.fasta"   2>/dev/null || echo 0)"
    echo "fungi_contigs=$(grep -c '^>' "${outdir}/contigs/fungi_contigs.fasta" 2>/dev/null || echo 0)"
    echo "other_contigs=$(grep -c '^>' "${outdir}/contigs/other_contigs.fasta" 2>/dev/null || echo 0)"
  } > "${outdir}/contigs/classify/tiara.counts.txt"
}

# ---------------------------
# 4) Run Prodigal (bacterial bucket)
# ---------------------------
run_prodigal_bac() {
  local sample_id="$1"
  local bac_fa="$2"
  local outdir="$3"

  mkdirp "${outdir}/bac/prodigal"
  if [[ ! -s "${bac_fa}" ]]; then
    msg "[${sample_id}] No bacterial contigs FASTA found (empty). Skipping Prodigal."
    return 0
  fi

  msg "[${sample_id}] Running Prodigal (meta mode) on bacterial contigs"
  prodigal \
    -i "${bac_fa}" \
    -p meta \
    -a "${outdir}/bac/prodigal/bac.proteins.faa" \
    -d "${outdir}/bac/prodigal/bac.cds.fna" \
    -f gff \
    -o "${outdir}/bac/prodigal/bac.genes.gff" \
    2> "${outdir}/bac/prodigal/prodigal.log"

  # Optional: bacterial small peptides (<=100 aa) from Prodigal output
  seqkit seq -M "${MAX_FUNGAL_PEPTIDE_AA}" "${outdir}/bac/prodigal/bac.proteins.faa" \
    > "${outdir}/bac/prodigal/bac.prodigal_le_${MAX_FUNGAL_PEPTIDE_AA}aa.faa"
}

# ---------------------------
# 5) Run SmORFinder (bacterial bucket)
# ---------------------------
run_smorfinder_bac() {
  local sample_id="$1"
  local bac_fa="$2"
  local outdir="$3"

  mkdirp "${outdir}/bac/smorfinder"
  if [[ ! -s "${bac_fa}" ]]; then
    msg "[${sample_id}] No bacterial contigs FASTA found (empty). Skipping SmORFinder."
    return 0
  fi

  msg "[${sample_id}] Running SmORFinder (meta) on bacterial contigs"
  # SmORFinder quickstart: smorf meta myMetagenome.fna :contentReference[oaicite:11]{index=11}
  # It also supports downloading needed data via `smorf` (no args) once installed. :contentReference[oaicite:12]{index=12}
  ( cd "${outdir}/bac/smorfinder" && smorf meta "../../contigs/bac_contigs.fasta" --threads "${CPUS}" ) \
    > "${outdir}/bac/smorfinder/smorfinder.stdout.log" \
    2> "${outdir}/bac/smorfinder/smorfinder.stderr.log" || true

  # We don't assume an exact output filename; we’ll later "find" peptide FASTAs for pooling.
}

# ---------------------------
# 6) Run funannotate (fungal bucket) + flag <=100 aa peptides
# ---------------------------
run_funannotate_fungi() {
  local sample_id="$1"
  local fungi_fa="$2"
  local outdir="$3"

  mkdirp "${outdir}/fungi"
  if [[ ! -s "${fungi_fa}" ]]; then
    msg "[${sample_id}] No fungal/euk contigs FASTA found (empty). Skipping funannotate."
    return 0
  fi

  msg "[${sample_id}] Running funannotate predict on fungal contigs (minimal run)"
  msg "[${sample_id}] NOTE: funannotate often needs databases set up (funannotate setup)."
  if [[ -n "${FUNANNOTATE_DB_DIR}" ]]; then
    msg "[${sample_id}] Using FUNANNOTATE_DB_DIR=${FUNANNOTATE_DB_DIR}"
    export FUNANNOTATE_DB="${FUNANNOTATE_DB_DIR}"
  fi

  # Minimal example: funannotate predict -i genome.fasta -o out --species "Name" --cpus N :contentReference[oaicite:13]{index=13}
  funannotate predict \
    -i "${fungi_fa}" \
    -o "${outdir}/fungi/funannotate_out" \
    --species "Fungus_sp_${sample_id}" \
    --cpus "${CPUS}" \
    > "${outdir}/fungi/funannotate.predict.stdout.log" \
    2> "${outdir}/fungi/funannotate.predict.stderr.log" || true

  # Find predicted protein FASTA (funannotate output naming can vary by version/config)
  msg "[${sample_id}] Locating funannotate predicted proteins FASTA"
  local prot_fa
  prot_fa="$(find "${outdir}/fungi/funannotate_out" -maxdepth 3 -type f \
              \( -iname "*protein*.fa" -o -iname "*protein*.faa" -o -iname "*proteins*.fa" -o -iname "*proteins*.faa" \) \
              | head -n 1 || true)"

  if [[ -z "${prot_fa}" ]]; then
    msg "[${sample_id}] WARNING: Could not find a proteins FASTA under funannotate output. Check logs."
    return 0
  fi

  msg "[${sample_id}] Flagging fungal peptides <= ${MAX_FUNGAL_PEPTIDE_AA} aa from: ${prot_fa}"
  seqkit seq -M "${MAX_FUNGAL_PEPTIDE_AA}" "${prot_fa}" \
    > "${outdir}/fungi/fungi_le_${MAX_FUNGAL_PEPTIDE_AA}aa.faa"
}

# ---------------------------
# 7) Pool peptides + de-replicate (non-redundant catalog)
# ---------------------------
create_nonredundant_catalog() {
  local sample_id="$1"
  local outdir="$2"

  mkdirp "${outdir}/catalog"
  local pooled="${outdir}/catalog/pooled.peptides.faa"

  msg "[${sample_id}] Pooling peptides (bacterial + fungal) into one FASTA"

  # Collect:
  #   - bacterial Prodigal <=100aa (baseline)
  #   - any SmORFinder-produced peptide/protein FASTAs (best-effort glob)
  #   - fungal <=100aa (from funannotate output)
  : > "${pooled}"

  if [[ -s "${outdir}/bac/prodigal/bac.prodigal_le_${MAX_FUNGAL_PEPTIDE_AA}aa.faa" ]]; then
    cat "${outdir}/bac/prodigal/bac.prodigal_le_${MAX_FUNGAL_PEPTIDE_AA}aa.faa" >> "${pooled}"
  fi

  # SmORFinder: best-effort find FASTA outputs inside its folder
  # (We do not hardcode filenames because versions may differ.)
  find "${outdir}/bac/smorfinder" -type f \( -iname "*.faa" -o -iname "*.fa" -o -iname "*.fasta" \) \
    | while read -r f; do
        # Skip the original contig fasta and other non-peptide files if present:
        case "$(basename "$f")" in
          *contigs*.fa* ) continue ;;
        esac
        # Heuristic: only append if it looks like protein FASTA (contains '>')
        if grep -q '^>' "$f"; then
          cat "$f" >> "${pooled}"
        fi
      done

  if [[ -s "${outdir}/fungi/fungi_le_${MAX_FUNGAL_PEPTIDE_AA}aa.faa" ]]; then
    cat "${outdir}/fungi/fungi_le_${MAX_FUNGAL_PEPTIDE_AA}aa.faa" >> "${pooled}"
  fi

  # If pooled is empty, skip.
  if [[ ! -s "${pooled}" ]]; then
    msg "[${sample_id}] No peptides found to pool. Skipping de-replication."
    return 0
  fi

  msg "[${sample_id}] De-replicating pooled peptides with MMseqs2 easy-cluster"
  # mmseqs easy-cluster input.fasta outPrefix tmpDir :contentReference[oaicite:14]{index=14}
  # Parameters are tunable; here is a reasonable starting point:
  #   --min-seq-id 0.95  (95% identity)
  #   -c 0.8             (80% coverage)
  #   --cov-mode 1       (coverage of target)
  mmseqs easy-cluster "${pooled}" "${outdir}/catalog/nr_catalog" "${outdir}/catalog/tmp" \
    --min-seq-id 0.95 -c 0.8 --cov-mode 1 \
    > "${outdir}/catalog/mmseqs.easycluster.stdout.log" \
    2> "${outdir}/catalog/mmseqs.easycluster.stderr.log" || true

  msg "[${sample_id}] Non-redundant outputs (representatives + clusters) are under: ${outdir}/catalog/nr_catalog*"
}

# ---------------------------
# Main per-sample runner
# ---------------------------
run_one_sample() {
  local sample_id="$1"
  local asm_fa="${ASSEMBLY_ROOT}/${sample_id}/assembly.fasta"
  [[ -s "${asm_fa}" ]] || die "Missing assembly for sample '${sample_id}': ${asm_fa}"

  local outdir="${RESULTS_DIR}/${sample_id}"
  mkdirp "${outdir}"

  flag_contigs_and_record_stats "${sample_id}" "${asm_fa}" "${outdir}"
  classify_contigs_tiara         "${sample_id}" "${asm_fa}" "${outdir}"
  run_prodigal_bac              "${sample_id}" "${outdir}/contigs/bac_contigs.fasta"   "${outdir}"
  run_smorfinder_bac            "${sample_id}" "${outdir}/contigs/bac_contigs.fasta"   "${outdir}"
  run_funannotate_fungi         "${sample_id}" "${outdir}/contigs/fungi_contigs.fasta" "${outdir}"
  create_nonredundant_catalog   "${sample_id}" "${outdir}"

  msg "[${sample_id}] DONE"
}

# ---------------------------
# CLI
# ---------------------------
usage() {
  cat <<EOF
Usage:
  workflow/run_smorfs_pipeline.sh --create-env
  workflow/run_smorfs_pipeline.sh --run --sample TS-0500 --cpus 8
  workflow/run_smorfs_pipeline.sh --run --samples-file samples.txt --cpus 16

Notes:
  - Expects MetaFlye assembly at: ${ASSEMBLY_ROOT}/<sample_id>/assembly.fasta
  - Appends to:
      ${RESULTS_DIR}/assembly_stats.tsv
      ${RESULTS_DIR}/contig_stats.tsv
  - Optional (recommended) funannotate DB:
      export FUNANNOTATE_DB_DIR=/path/to/funannotate_db
      funannotate setup -d \$FUNANNOTATE_DB_DIR

Args:
  --create-env
  --run
  --sample <id>
  --samples-file <file>     (one sample_id per line)
  --cpus <int>              (default: ${CPUS})
  --min-contig-bp <int>     (default: ${MIN_CONTIG_BP})
  --max-fungal-aa <int>     (default: ${MAX_FUNGAL_PEPTIDE_AA})
EOF
}

MODE=""
SAMPLE=""
SAMPLES_FILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --create-env) MODE="create-env"; shift ;;
    --run)        MODE="run"; shift ;;
    --sample)     SAMPLE="${2:-}"; shift 2 ;;
    --samples-file) SAMPLES_FILE="${2:-}"; shift 2 ;;
    --cpus)       CPUS="${2:-}"; shift 2 ;;
    --min-contig-bp) MIN_CONTIG_BP="${2:-}"; shift 2 ;;
    --max-fungal-aa) MAX_FUNGAL_PEPTIDE_AA="${2:-}"; shift 2 ;;
    -h|--help)    usage; exit 0 ;;
    *) die "Unknown argument: $1 (use --help)" ;;
  esac
done

if [[ -z "${MODE}" ]]; then
  usage
  exit 1
fi

if [[ "${MODE}" == "create-env" ]]; then
  create_env
  exit 0
fi

# MODE=run
ensure_env_once
init_tables

need_cmd seqkit
need_cmd prodigal
need_cmd mmseqs
need_cmd tiara
need_cmd smorf || msg "WARNING: 'smorf' not found. Did you install SmORFinder in the env?"
need_cmd funannotate || msg "WARNING: 'funannotate' not found. Did you install funannotate in the env?"

if [[ -n "${SAMPLE}" ]]; then
  run_one_sample "${SAMPLE}"
elif [[ -n "${SAMPLES_FILE}" ]]; then
  [[ -f "${SAMPLES_FILE}" ]] || die "Samples file not found: ${SAMPLES_FILE}"
  while read -r sid; do
    [[ -z "${sid}" ]] && continue
    run_one_sample "${sid}"
  done < "${SAMPLES_FILE}"
else
  die "Provide --sample <id> or --samples-file <file>"
fi
