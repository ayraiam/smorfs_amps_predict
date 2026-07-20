#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# Workflow controls
# ============================================================

# 1 = run the initial smORF filtering and GFF annotation
# 0 = skip this step
RUN_INITIAL_FILTER="${RUN_INITIAL_FILTER:-1}"

# 1 = run the Prodigal completeness filtering and ranking
# 0 = skip this step
RUN_PRODIGAL_RANKING="${RUN_PRODIGAL_RANKING:-1}"


# ============================================================
# Parse command-line flags
# ============================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --rank-only)
            RUN_INITIAL_FILTER=0
            RUN_PRODIGAL_RANKING=1
            shift
            ;;

        --initial-only)
            RUN_INITIAL_FILTER=1
            RUN_PRODIGAL_RANKING=0
            shift
            ;;

        --help|-h)
            cat <<'EOF'
Usage:

  bash workflow/run_filter_and_annotate_smorfs.sh [option]

Options:

  --rank-only
      Skip the initial filtering/GFF annotation step and run only
      the Prodigal completeness filtering and ranking step.

  --initial-only
      Run only the initial filtering/GFF annotation step.

  --help, -h
      Show this help message.

With no option, both steps are run.
EOF
            exit 0
            ;;

        *)
            echo "ERROR: Unknown option: $1" >&2
            echo "Use --help to see the available options." >&2
            exit 1
            ;;
    esac
done


# ============================================================
# Scripts
# ============================================================

FILTER_SCRIPT="${FILTER_SCRIPT:-workflow/filter_and_annotate_smorfs.R}"

RANK_SCRIPT="${RANK_SCRIPT:-workflow/rank_smorfs_by_prodigal.R}"


# ============================================================
# Input/output paths
# ============================================================

INPUT_TSV="${INPUT_TSV:-results/smorfs/UNIFORME_GLOBAL/catalog/predicted_smorfs.with_macrel.refined_bacs.with_clusters.tsv}"

INPUT_GFF="${INPUT_GFF:-results/smorfs/UNIFORME_GLOBAL/bac/prodigal/bac.genes.gff}"

SECONDARY_GFF="${SECONDARY_GFF:-results/smorfs/UNIFORME_GLOBAL/bac/smorfinder/smorf_output/smorf_output.gff}"

INITIAL_FILTERED_TSV="${INITIAL_FILTERED_TSV:-results/smorfs/UNIFORME_GLOBAL/catalog/predicted_smorfs.initial_filtered.with_prodigal_gff.tsv}"

RANKED_OUTPUT_TSV="${RANKED_OUTPUT_TSV:-results/smorfs/UNIFORME_GLOBAL/catalog/predicted_smorfs.complete_ranked_by_prodigal.tsv}"


# ============================================================
# Initial filtering parameters
# ============================================================

MATCHING_COLUMN="${MATCHING_COLUMN:-feature_id}"

MAX_AA_LEN="${MAX_AA_LEN:-100}"

MIN_AMP_PROB="${MIN_AMP_PROB:-0.5}"


# ============================================================
# Logging
# ============================================================

LOGDIR="${LOGDIR:-logs}"
mkdir -p "$LOGDIR"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

OUT_LOG="${OUT_LOG:-${LOGDIR}/filter_and_rank_${TIMESTAMP}.out}"
ERR_LOG="${ERR_LOG:-${LOGDIR}/filter_and_rank_${TIMESTAMP}.err}"
COMMAND_LOG="${COMMAND_LOG:-${LOGDIR}/filter_and_rank_run_command_${TIMESTAMP}.txt}"

{
    echo "Command executed:"
    printf '%q ' "$0" "$@"
    echo
    echo

    echo "Working directory:"
    pwd
    echo

    echo "Date:"
    date
    echo

    echo "Hostname:"
    hostname
    echo

    echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-NA}"
    echo

    echo "Workflow controls:"
    echo "RUN_INITIAL_FILTER=$RUN_INITIAL_FILTER"
    echo "RUN_PRODIGAL_RANKING=$RUN_PRODIGAL_RANKING"
    echo

    echo "Scripts:"
    echo "FILTER_SCRIPT=$FILTER_SCRIPT"
    echo "RANK_SCRIPT=$RANK_SCRIPT"
    echo

    echo "Input/output paths:"
    echo "INPUT_TSV=$INPUT_TSV"
    echo "INPUT_GFF=$INPUT_GFF"
    echo "SECONDARY_GFF=$SECONDARY_GFF"
    echo "INITIAL_FILTERED_TSV=$INITIAL_FILTERED_TSV"
    echo "RANKED_OUTPUT_TSV=$RANKED_OUTPUT_TSV"
    echo

    echo "Filtering parameters:"
    echo "MATCHING_COLUMN=$MATCHING_COLUMN"
    echo "MAX_AA_LEN=$MAX_AA_LEN"
    echo "MIN_AMP_PROB=$MIN_AMP_PROB"
} > "$COMMAND_LOG"

exec >"$OUT_LOG" 2>"$ERR_LOG"

echo "Command log:        $COMMAND_LOG"
echo "Standard output:    $OUT_LOG"
echo "Standard error:     $ERR_LOG"


# ============================================================
# Initial filtering and GFF annotation
# ============================================================

if [[ "$RUN_INITIAL_FILTER" == "1" ]]; then

    echo
    echo "============================================================"
    echo "Running initial smORF filtering and GFF annotation"
    echo "============================================================"
    echo "Input TSV:          $INPUT_TSV"
    echo "Primary GFF:        $INPUT_GFF"
    echo "Secondary GFF:      $SECONDARY_GFF"
    echo "Output TSV:         $INITIAL_FILTERED_TSV"
    echo "Matching column:    $MATCHING_COLUMN"
    echo "Maximum aa_len:     $MAX_AA_LEN"
    echo "Minimum amp_prob:   $MIN_AMP_PROB"
    echo "Required flag_edge: 0"
    echo "============================================================"

    if [[ ! -f "$FILTER_SCRIPT" ]]; then
        echo "ERROR: Filtering R script not found: $FILTER_SCRIPT" >&2
        exit 1
    fi

    if [[ ! -f "$INPUT_TSV" ]]; then
        echo "ERROR: Input TSV not found: $INPUT_TSV" >&2
        exit 1
    fi

    if [[ ! -f "$INPUT_GFF" ]]; then
        echo "ERROR: Primary GFF not found: $INPUT_GFF" >&2
        exit 1
    fi

    if [[ ! -f "$SECONDARY_GFF" ]]; then
        echo "ERROR: Secondary GFF not found: $SECONDARY_GFF" >&2
        exit 1
    fi

    mkdir -p "$(dirname "$INITIAL_FILTERED_TSV")"

    Rscript "$FILTER_SCRIPT" \
        "$INPUT_TSV" \
        "$INPUT_GFF" \
        "$SECONDARY_GFF" \
        "$INITIAL_FILTERED_TSV" \
        "$MATCHING_COLUMN" \
        "$MAX_AA_LEN" \
        "$MIN_AMP_PROB"

    echo
    echo "Initial filtering completed:"
    echo "  $INITIAL_FILTERED_TSV"

else

    echo
    echo "Skipping initial smORF filtering and GFF annotation."

fi


# ============================================================
# Prodigal completeness filtering and ranking
# ============================================================

if [[ "$RUN_PRODIGAL_RANKING" == "1" ]]; then

    echo
    echo "============================================================"
    echo "Running Prodigal completeness filtering and ranking"
    echo "============================================================"
    echo "Input TSV:  $INITIAL_FILTERED_TSV"
    echo "Output TSV: $RANKED_OUTPUT_TSV"
    echo
    echo "Filtering:"
    echo "  gff_attr_partial == 00"
    echo
    echo "Ranking:"
    echo "  1. gff_attr_conf   decreasing"
    echo "  2. gff_attr_score  decreasing"
    echo "  3. gff_attr_cscore decreasing"
    echo "============================================================"

    if [[ ! -f "$RANK_SCRIPT" ]]; then
        echo "ERROR: Ranking R script not found: $RANK_SCRIPT" >&2
        exit 1
    fi

    if [[ ! -s "$INITIAL_FILTERED_TSV" ]]; then
        echo "ERROR: Annotated TSV not found or empty:" >&2
        echo "  $INITIAL_FILTERED_TSV" >&2
        echo >&2
        echo "Run the initial step first or verify the path." >&2
        exit 1
    fi

    mkdir -p "$(dirname "$RANKED_OUTPUT_TSV")"

    Rscript "$RANK_SCRIPT" \
        "$INITIAL_FILTERED_TSV" \
        "$RANKED_OUTPUT_TSV"

    echo
    echo "Prodigal ranking completed:"
    echo "  $RANKED_OUTPUT_TSV"

else

    echo
    echo "Skipping Prodigal completeness filtering and ranking."

fi


echo
echo "============================================================"
echo "Workflow completed successfully."
echo "============================================================"
