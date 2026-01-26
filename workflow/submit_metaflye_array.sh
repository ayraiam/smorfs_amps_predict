#!/usr/bin/env bash
# ==========================================================
# Script: workflow/submit_metaflye_array.sh
# Purpose: Batch-safe MetaFlye co-assembly submission.
#   - Scans FASTQs present in FASTQ_DIR (default: data/)
#   - For each FASTQ present, requires a SampleID mapping in METADATA_MAP
#   - Builds a SampleID->FASTQ map for ONLY present FASTQs
#   - Submits Slurm array: 1 task per SampleID
#
# Optional:
#   --strict-metadata
#     Also require that every FASTQ listed in METADATA_MAP exists in FASTQ_DIR.
# ==========================================================
set -euo pipefail

PARTITION="${PARTITION:-short}"
TIME="${TIME:-12:00:00}"
CPUS="${CPUS:-16}"
MEM="${MEM:-64G}"
WDIR="${WDIR:-$PWD}"

FASTQ_DIR="${FASTQ_DIR:-data}"
RESULTS_DIR="${RESULTS_DIR:-results}"

# User-provided mapping file: SampleID <-> FASTQ filename/path
METADATA_MAP="${METADATA_MAP:-metadata/metagenome_files.txt}"

# Optional knobs exported to tasks
READ_MODE="${READ_MODE:-nano-raw}"     # nano-raw recommended
MIN_OVERLAP="${MIN_OVERLAP:-}"         # empty => auto
GENOME_SIZE="${GENOME_SIZE:-}"         # usually empty for metagenomes

mkdir -p logs metadata

TS="$(date +%Y%m%d_%H%M%S)"
SAMPLE_LIST="metadata/metaflye_sampleids_${TS}.list"
MAP_PRESENT="metadata/metaflye_sample_fastqs_${TS}.tsv"

STRICT_METADATA=0

die() { echo "ERROR: $*" >&2; exit 1; }

usage() {
  cat <<EOF
Usage:
  bash workflow/submit_metaflye_array.sh [--strict-metadata]

Batch-safe default behavior:
  1) Scan FASTQs present in FASTQ_DIR (${FASTQ_DIR})
  2) For each FASTQ present, look up its SampleID in METADATA_MAP (${METADATA_MAP})
  3) If a FASTQ in FASTQ_DIR is missing from metadata -> ERROR (likely typo/missing row)
  4) Ignore metadata rows for FASTQs not present yet (supports running in batches)

Options:
  --strict-metadata
    Also require that every FASTQ listed in METADATA_MAP exists in FASTQ_DIR.

EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strict-metadata) STRICT_METADATA=1; shift 1 ;;
    -h|--help) usage ;;
    *) die "Unknown argument: $1 (try --help)" ;;
  esac
done

detect_delim() {
  local f="$1"
  if head -n 1 "$f" | grep -q $'\t'; then
    echo $'\t'
  elif head -n 1 "$f" | grep -q ','; then
    echo ','
  else
    echo ' '
  fi
}

normalize_map_to_tsv() {
  # Normalized output: SampleID<TAB>FASTQ_FIELD (as provided)
  local in="$1"
  local out="$2"
  [[ -f "$in" ]] || die "Metadata map not found: $in"

  local delim header
  delim="$(detect_delim "$in")"
  header="$(head -n 1 "$in")"

  local SAMPLE_COL_SEL="${SAMPLE_COL:-}"
  local FASTQ_COL_SEL="${FASTQ_COL:-}"

  # Helper: resolve column name -> index (1-based)
  colname_to_idx() {
    local name="$1"
    awk -v FS="$delim" -v target="$(echo "$name" | tr '[:upper:]' '[:lower:]')" '
      NR==1{
        for(i=1;i<=NF;i++){
          h=tolower($i)
          gsub(/^[[:space:]]+|[[:space:]]+$/,"",h)
          if(h==target){ print i; exit }
        }
      }
    ' "$in"
  }

  # If user specified columns, use them
  if [[ -n "$SAMPLE_COL_SEL" && -n "$FASTQ_COL_SEL" ]]; then
    local sidx fidx

    if [[ "$SAMPLE_COL_SEL" =~ ^[0-9]+$ ]]; then
      sidx="$SAMPLE_COL_SEL"
    else
      sidx="$(colname_to_idx "$SAMPLE_COL_SEL")"
    fi

    if [[ "$FASTQ_COL_SEL" =~ ^[0-9]+$ ]]; then
      fidx="$FASTQ_COL_SEL"
    else
      fidx="$(colname_to_idx "$FASTQ_COL_SEL")"
    fi

    [[ -n "${sidx}" ]] || die "Could not resolve SAMPLE_COL=${SAMPLE_COL_SEL} to a column index"
    [[ -n "${fidx}" ]] || die "Could not resolve FASTQ_COL=${FASTQ_COL_SEL} to a column index"

    awk -v FS="$delim" -v OFS="\t" -v s="$sidx" -v f="$fidx" '
      NR==1{ next }  # skip header
      {
        sid=$s; fq=$f
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",sid)
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",fq)
        if(sid=="" || fq=="") next
        print sid, fq
      }
    ' "$in" > "$out"

    [[ -s "$out" ]] || die "Normalized map is empty after using SAMPLE_COL/FASTQ_COL. Check file format."
    return 0
  fi

  # Otherwise: auto-detect (your previous logic)
  if echo "$header" | grep -qiE 'sample' && echo "$header" | grep -qiE 'fastq|filename|file'; then
    awk -v FS="$delim" -v OFS="\t" '
      NR==1{
        for(i=1;i<=NF;i++){
          h=tolower($i)
          gsub(/^[[:space:]]+|[[:space:]]+$/,"",h)
          if(h ~ /^sampleid$/ || h ~ /sample/){ s=i }
          if(h ~ /fastq/ || h ~ /filename/ || h ~ /^file$/){ f=i }
        }
        if(!s || !f){
          print "ERROR: Could not find SampleID/FASTQ columns in header" > "/dev/stderr"
          exit 2
        }
        next
      }
      {
        sid=$s; fq=$f
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",sid)
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",fq)
        if(sid=="" || fq=="") next
        print sid, fq
      }
    ' "$in" > "$out"
  else
    awk -v FS="$delim" -v OFS="\t" '
      {
        sid=$1; fq=$2
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",sid)
        gsub(/^[[:space:]]+|[[:space:]]+$/,"",fq)
        if(sid=="" || fq=="") next
        print sid, fq
      }
    ' "$in" > "$out"
  fi

  [[ -s "$out" ]] || die "Normalized map is empty (check $in)"
}

list_fastqs_present_abs() {
  # Absolute paths for FASTQs currently present in FASTQ_DIR
  local -a files=()
  shopt -s nullglob
  files+=( "${FASTQ_DIR}"/*.fastq "${FASTQ_DIR}"/*.fastq.gz "${FASTQ_DIR}"/*.fq "${FASTQ_DIR}"/*.fq.gz )
  shopt -u nullglob

  [[ "${#files[@]}" -gt 0 ]] || die "No FASTQ files found in ${FASTQ_DIR}/"

  for f in "${files[@]}"; do
    if [[ "$f" = /* ]]; then
      echo "$f"
    else
      echo "${WDIR}/${f}"
    fi
  done | sort -V
}

build_present_only_map() {
  # Output: SampleID<TAB>ABS_FASTQ for FASTQs present now
  # Enforce:
  #   - Every FASTQ present must be in metadata (else error)
  # Optional strict:
  #   - Every metadata FASTQ must exist in FASTQ_DIR (else error)
  local norm_map="$1"
  local out="$2"

  : > "$out"

  # LUT: basename -> SampleID, plus keep original fq field for debug
  local lut="metadata/.tmp_lut_${TS}.tsv"
  awk -F $'\t' -v OFS="\t" '
    {
      sid=$1; fq=$2
      gsub(/^[[:space:]]+|[[:space:]]+$/,"",sid)
      gsub(/^[[:space:]]+|[[:space:]]+$/,"",fq)
      if(sid=="" || fq=="") next
      bn=fq
      sub(/^.*\//,"",bn)
      print bn, sid, fq
    }
  ' "$norm_map" > "$lut"

  # Detect duplicate basenames in metadata => ambiguous
  local dups
  dups="$(cut -f1 "$lut" | sort | uniq -d || true)"
  if [[ -n "$dups" ]]; then
    echo "ERROR: Duplicate FASTQ basenames in ${METADATA_MAP} (ambiguous mapping):" >&2
    echo "$dups" | sed 's/^/  - /' >&2
    echo "Fix: make FASTQ basenames unique, or use unique names in metadata." >&2
    exit 1
  fi

  local present_count=0 mapped_count=0
  local -a missing=()

  while IFS= read -r absf; do
    present_count=$((present_count+1))
    bn="$(basename "$absf")"

    sid="$(awk -F $'\t' -v b="$bn" '$1==b {print $2; exit}' "$lut")"
    if [[ -z "${sid}" ]]; then
      missing+=( "$bn" )
      continue
    fi

    mapped_count=$((mapped_count+1))
    echo -e "${sid}\t${absf}" >> "$out"
  done < <(list_fastqs_present_abs)

  # If any FASTQ present has no metadata mapping => hard fail (typo/missing row)
  if [[ "${#missing[@]}" -gt 0 ]]; then
    echo "ERROR: FASTQ file(s) present in ${FASTQ_DIR}/ but missing from ${METADATA_MAP}:" >&2
    for m in "${missing[@]}"; do echo "  - ${m}" >&2; done
    echo >&2
    echo "Likely causes: typo in filename, missing metadata row, wrong delimiter, or wrong column." >&2
    exit 1
  fi

  [[ -s "$out" ]] || die "No FASTQs mapped after scanning ${FASTQ_DIR}/"

  # Strict mode: metadata should not contain FASTQs that aren't present in FASTQ_DIR
  if [[ "${STRICT_METADATA}" -eq 1 ]]; then
    local present_basenames="metadata/.tmp_present_basenames_${TS}.txt"
    local meta_basenames="metadata/.tmp_meta_basenames_${TS}.txt"

    list_fastqs_present_abs | awk -F/ '{print $NF}' | sort -u > "$present_basenames"
    cut -f1 "$lut" | sort -u > "$meta_basenames"

    local not_present
    not_present="$(comm -23 "$meta_basenames" "$present_basenames" || true)"
    if [[ -n "$not_present" ]]; then
      echo "ERROR (--strict-metadata): metadata lists FASTQ(s) not present in ${FASTQ_DIR}/:" >&2
      echo "$not_present" | sed 's/^/  - /' >&2
      exit 1
    fi
  fi

  rm -f "$lut" "metadata/.tmp_present_basenames_${TS}.txt" "metadata/.tmp_meta_basenames_${TS}.txt" 2>/dev/null || true
}

build_sample_list() {
  local map="$1"
  local out="$2"
  cut -f1 "$map" | sort -u > "$out"
  [[ -s "$out" ]] || die "Sample list is empty"
}

submit_array() {
  local sample_list="$1"
  local map_file="$2"

  local n fastq_n
  n="$(wc -l < "$sample_list" | tr -d ' ')"
  fastq_n="$(wc -l < "$map_file" | tr -d ' ')"
  [[ "$n" -ge 1 ]] || die "No SampleIDs to submit"

  echo "MetaFlye submission summary: FASTQs_present=${fastq_n} | SampleIDs=${n} | array_tasks=${n} | strict=${STRICT_METADATA}"
  echo "Sample list: ${sample_list}"
  echo "Sample->FASTQ map (present only): ${map_file}"

  local OUT_LOG="logs/metaflye_array_${TS}.%A_%a.out"
  local ERR_LOG="logs/metaflye_array_${TS}.%A_%a.err"

  local jid
  jid="$(
    sbatch \
      --partition="${PARTITION}" \
      --time="${TIME}" \
      --cpus-per-task="${CPUS}" \
      --mem="${MEM}" \
      --chdir="${WDIR}" \
      --job-name="metaflye" \
      --array="1-${n}" \
      --output="${OUT_LOG}" \
      --error="${ERR_LOG}" \
      --export=ALL,THREADS="${CPUS}",RESULTS_DIR="${RESULTS_DIR}",METAFlyE_SAMPLE_LIST="${sample_list}",METAFlyE_SAMPLE_MAP="${map_file}",READ_MODE="${READ_MODE}",MIN_OVERLAP="${MIN_OVERLAP}",GENOME_SIZE="${GENOME_SIZE}" \
      workflow/metaflye_array_task.sh \
      | awk '{print $4}'
  )"

  echo "Submitted MetaFlye array job: ${jid}"
  echo "Per-task logs:"
  echo "  ${OUT_LOG}"
  echo "  ${ERR_LOG}"
}

# --------------------------- Main ---------------------------
TMP_NORM="metadata/.tmp_norm_${TS}.tsv"

normalize_map_to_tsv "${METADATA_MAP}" "${TMP_NORM}"
build_present_only_map "${TMP_NORM}" "${MAP_PRESENT}"
rm -f "${TMP_NORM}" || true

build_sample_list "${MAP_PRESENT}" "${SAMPLE_LIST}"
submit_array "${SAMPLE_LIST}" "${MAP_PRESENT}"
