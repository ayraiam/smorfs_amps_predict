#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Filter annotated smORF candidates to complete CDSs and rank
# them using Prodigal metrics.
#
# Input:
#   TSV produced by filter_and_annotate_smorfs.R
#
# Required columns:
#   gff_attr_partial
#   gff_attr_conf:
# Prodigal confidence score  - An overall confidence estimate assigned by Prodigal to the predicted CDS, 
# integrating multiple sequence features to reflect the reliability of the gene prediction. Higher values 
# indicate greater confidence that the predicted CDS represents a true protein-coding gene
#   gff_attr_score
# Overall Prodigal gene score - The total score assigned by Prodigal to the predicted CDS, representing the 
# combined evidence supporting the gene prediction. This score is approximately the sum of the coding score 
# (cscore) and the start-site score (sscore), with higher values indicating stronger overall support for the 
# predicted gene
#   gff_attr_cscore
# Coding potential score - The component of the Prodigal gene score that evaluates the coding potential of the 
# predicted CDS independently of translation initiation signals. Higher values indicate stronger evidence that
# the sequence represents a genuine protein-coding region.
#
# Filtering:
#   Keep only complete CDSs:
#     gff_attr_partial == "00"
#
# Ranking:
#   1. gff_attr_conf   decreasing. 
#   2. gff_attr_score  decreasing, as tie breaker. 
#   3. gff_attr_cscore decreasing, as tie breaker.  
#
# Output:
#   All original input columns, followed by:
#     prodigal_rank
# ============================================================


# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop(
    paste0(
      "\nUsage:\n",
      "  Rscript rank_smorfs_by_prodigal.R ",
      "<input_annotated.tsv> <output_ranked.tsv>\n\n",
      "Example:\n",
      "  Rscript rank_smorfs_by_prodigal.R ",
      "results/smorfs/UNIFORME_GLOBAL/catalog/",
      "predicted_smorfs.initial_filtered.with_prodigal_gff.tsv ",
      "results/smorfs/UNIFORME_GLOBAL/catalog/",
      "predicted_smorfs.complete_ranked_by_prodigal.tsv\n"
    ),
    call. = FALSE
  )
}

input_tsv  <- args[1]
output_tsv <- args[2]


# ------------------------------------------------------------
# Validate paths
# ------------------------------------------------------------

if (!file.exists(input_tsv)) {
  stop(
    "Input annotated TSV does not exist: ",
    input_tsv,
    call. = FALSE
  )
}

output_dir <- dirname(output_tsv)

if (!dir.exists(output_dir)) {
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# Read input
# ------------------------------------------------------------

message("Reading annotated TSV:")
message("  ", input_tsv)

dt <- fread(
  input_tsv,
  sep = "\t",
  header = TRUE,
  quote = "",
  na.strings = c("", "NA"),
  data.table = TRUE,
  showProgress = TRUE
)

original_column_order <- names(dt)
n_input <- nrow(dt)

message("  Input rows: ", n_input)


# ------------------------------------------------------------
# Validate required columns
# ------------------------------------------------------------

required_columns <- c(
  "gff_attr_partial",
  "gff_attr_conf",
  "gff_attr_score",
  "gff_attr_cscore"
)

missing_columns <- setdiff(
  required_columns,
  names(dt)
)

if (length(missing_columns) > 0L) {
  stop(
    "The following required columns are missing: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Preserve original row order for deterministic tie handling
# ------------------------------------------------------------

dt[, original_input_order__ := .I]


# ------------------------------------------------------------
# Convert Prodigal metrics to numeric
# ------------------------------------------------------------

dt[, prodigal_conf__ := suppressWarnings(
  as.numeric(gff_attr_conf)
)]

dt[, prodigal_score__ := suppressWarnings(
  as.numeric(gff_attr_score)
)]

dt[, prodigal_cscore__ := suppressWarnings(
  as.numeric(gff_attr_cscore)
)]


# ------------------------------------------------------------
# Keep only complete CDSs
# ------------------------------------------------------------

complete_dt <- dt[
  !is.na(gff_attr_partial) &
    as.character(gff_attr_partial) %in% c("00", "0")
]

n_complete <- nrow(complete_dt)
n_partial_removed <- n_input - n_complete


# ------------------------------------------------------------
# Remove complete CDSs lacking ranking metrics
# ------------------------------------------------------------

valid_rank_dt <- complete_dt[
  !is.na(prodigal_conf__) &
    !is.na(prodigal_score__) &
    !is.na(prodigal_cscore__)
]

n_valid <- nrow(valid_rank_dt)
n_missing_metrics <- n_complete - n_valid


# ------------------------------------------------------------
# Print filtering summary
# ------------------------------------------------------------

message("")
message("Prodigal ranking summary:")
message("  Input annotated rows:             ", n_input)
message("  Complete CDS rows (partial == 0/00): ", n_complete)
message("  Partial CDS rows removed:          ", n_partial_removed)
message("  Complete rows missing metrics:     ", n_missing_metrics)
message("  Final rows available for ranking:  ", n_valid)


# ------------------------------------------------------------
# Handle empty output
# ------------------------------------------------------------

if (n_valid == 0L) {
  warning(
    "No complete CDSs with valid conf, score, and cscore ",
    "were available for ranking."
  )
  
  valid_rank_dt[, prodigal_rank := integer()]
  
  empty_result <- valid_rank_dt[
    ,
    c(
      original_column_order,
      "prodigal_rank"
    ),
    with = FALSE
  ]
  
  fwrite(
    empty_result,
    output_tsv,
    sep = "\t",
    quote = FALSE,
    na = ""
  )
  
  message("")
  message("Output written:")
  message("  ", output_tsv)
  
  quit(
    save = "no",
    status = 0
  )
}


# ------------------------------------------------------------
# Rank candidates
# ------------------------------------------------------------

setorder(
  valid_rank_dt,
  -prodigal_conf__,
  -prodigal_score__,
  -prodigal_cscore__,
  original_input_order__
)

valid_rank_dt[, prodigal_rank := .I]


# ------------------------------------------------------------
# Prepare final output
# ------------------------------------------------------------

final_output_columns <- c(
  original_column_order,
  "prodigal_rank"
)

result <- valid_rank_dt[
  ,
  final_output_columns,
  with = FALSE
]


# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

fwrite(
  result,
  output_tsv,
  sep = "\t",
  quote = FALSE,
  na = ""
)

message("")
message("Completed successfully.")
message("Output file:")
message("  ", output_tsv)
message("")
message("Output dimensions:")
message("  Rows:    ", nrow(result))
message("  Columns: ", ncol(result))
message("")
message("Top-ranked candidates:")

top_n <- min(10L, nrow(result))

for (i in seq_len(top_n)) {
  message(
    "  Rank ", result$prodigal_rank[i],
    ": conf=", result$gff_attr_conf[i],
    "; score=", result$gff_attr_score[i],
    "; cscore=", result$gff_attr_cscore[i]
  )
}
