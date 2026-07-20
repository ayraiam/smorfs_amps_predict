#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# Filter a smORF/Macrel TSV and append matching Prodigal GFF
# information.
#
# Filtering:
#   #   aa_len <= maximum length
#   amp_prob > minimum probability
#   flag_edge == 0
#
# Matching:
#   TSV feature_id: contig_1_9
#
#   GFF:
#     column 1: contig_1
#     column 9: ID=1_9;...
#
#   Constructed GFF key:
#     contig_1 + "_" + final component of ID
#     = contig_1_9
#
# Output:
#   All original TSV columns, in their original order,
#   followed by:
#     - the eight standard GFF columns
#     - the complete raw ninth GFF column
#     - one separate column for every GFF attribute
# ============================================================


# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop(
    paste0(
      "\nUsage:\n",
      "  Rscript filter_and_annotate_smorfs.R ",
      "<input.tsv> <primary.gff> <secondary.gff> <output.tsv> ",
      "[matching_column] [max_aa_len] [min_amp_prob]\n\n",
      "Example:\n",
      "  Rscript filter_and_annotate_smorfs.R ",
      "predicted_smorfs.with_macrel.refined_bacs.with_clusters.tsv ",
      "bac.genes.gff ",
      "smorf_output.gff ",
      "filtered_smorfs.with_gff.tsv ",
      "feature_id 100 0.5\n"
    ),
    call. = FALSE
  )
}

input_tsv       <- args[1]
input_gff       <- args[2]
secondary_gff   <- args[3]
output_tsv      <- args[4]
matching_col    <- if (length(args) >= 5) args[5] else "feature_id"
max_aa_len      <- if (length(args) >= 6) as.numeric(args[6]) else 100
min_amp_prob    <- if (length(args) >= 7) as.numeric(args[7]) else 0.5


# ------------------------------------------------------------
# Validate input files
# ------------------------------------------------------------

if (!file.exists(input_tsv)) {
  stop("Input TSV does not exist: ", input_tsv, call. = FALSE)
}

if (!file.exists(input_gff)) {
  stop("Input GFF does not exist: ", input_gff, call. = FALSE)
}

if (!file.exists(secondary_gff)) {
  stop(
    "Secondary GFF does not exist: ",
    secondary_gff,
    call. = FALSE
  )
}

if (is.na(max_aa_len)) {
  stop("max_aa_len must be numeric.", call. = FALSE)
}

if (is.na(min_amp_prob)) {
  stop("min_amp_prob must be numeric.", call. = FALSE)
}

output_dir <- dirname(output_tsv)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

# Extract a particular key from a semicolon-delimited GFF
# attributes string.
extract_gff_attribute <- function(attribute_vector, key) {
  pattern <- paste0(
    "(?:^|;)",
    key,
    "=([^;]*)"
  )
  
  result <- rep(NA_character_, length(attribute_vector))
  
  matched <- regexec(
    pattern,
    attribute_vector,
    perl = TRUE
  )
  
  extracted <- regmatches(attribute_vector, matched)
  
  has_match <- lengths(extracted) >= 2L
  
  result[has_match] <- vapply(
    extracted[has_match],
    function(x) x[2],
    character(1)
  )
  
  result
}


# Return all attribute names occurring in GFF column 9.
get_all_attribute_names <- function(attribute_vector) {
  attribute_vector <- attribute_vector[
    !is.na(attribute_vector) &
      nzchar(attribute_vector) &
      attribute_vector != "."
  ]
  
  if (length(attribute_vector) == 0L) {
    return(character())
  }
  
  attribute_parts <- unlist(
    strsplit(attribute_vector, ";", fixed = TRUE),
    use.names = FALSE
  )
  
  attribute_parts <- attribute_parts[nzchar(attribute_parts)]
  
  attribute_names <- sub("=.*$", "", attribute_parts)
  
  attribute_names <- trimws(attribute_names)
  attribute_names <- attribute_names[nzchar(attribute_names)]
  
  sort(unique(attribute_names))
}

read_and_parse_gff <- function(
    gff_path,
    annotation_source,
    match_mode = c("primary", "exact_id")
) {
  match_mode <- match.arg(match_mode)
  
  message("")
  message("Reading GFF:")
  message("  ", gff_path)
  
  awk_command <- sprintf(
    "awk 'BEGIN {FS=\"\\t\"} $0 !~ /^#/ && NF >= 9 {print}' %s",
    shQuote(gff_path)
  )
  
  gff_dt <- fread(
    cmd = awk_command,
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    na.strings = character(),
    data.table = TRUE,
    showProgress = TRUE
  )
  
  if (ncol(gff_dt) < 9L) {
    stop(
      "The parsed GFF has fewer than nine columns: ",
      gff_path,
      call. = FALSE
    )
  }
  
  if (ncol(gff_dt) > 9L) {
    warning(
      "The GFF produced more than nine parsed columns: ",
      gff_path,
      ". Only the first nine columns will be used."
    )
    
    gff_dt <- gff_dt[, 1:9]
  }
  
  setnames(
    gff_dt,
    c(
      "gff_seqid",
      "gff_source",
      "gff_type",
      "gff_start",
      "gff_end",
      "gff_score",
      "gff_strand",
      "gff_phase",
      "gff_attributes_raw"
    )
  )
  
  message("  Parsed GFF feature rows:     ", nrow(gff_dt))
  message("  Parsing GFF attributes...")
  
  attribute_names <- get_all_attribute_names(
    gff_dt$gff_attributes_raw
  )
  
  message(
    "  Distinct attribute fields:   ",
    length(attribute_names)
  )
  
  if (length(attribute_names) > 0L) {
    for (attribute_name in attribute_names) {
      output_attribute_name <- paste0(
        "gff_attr_",
        attribute_name
      )
      
      gff_dt[
        ,
        (output_attribute_name) :=
          extract_gff_attribute(
            gff_attributes_raw,
            attribute_name
          )
      ]
    }
  }
  
  if (!"gff_attr_ID" %in% names(gff_dt)) {
    stop(
      "No ID attribute was found in GFF column 9: ",
      gff_path,
      call. = FALSE
    )
  }
  
  if (match_mode == "primary") {
    gff_dt[
      ,
      gff_id_final_component__ :=
        sub("^.*_", "", gff_attr_ID)
    ]
    
    gff_dt[
      ,
      gff_match_key__ :=
        paste0(
          gff_seqid,
          "_",
          gff_id_final_component__
        )
    ]
  } else {
    # Secondary GFF IDs are matched directly against feature_id.
    gff_dt[, gff_id_final_component__ := NA_character_]
    gff_dt[, gff_match_key__ := as.character(gff_attr_ID)]
  }
  
  gff_dt[, gff_annotation_source := annotation_source]
  
  gff_dt
}


# ------------------------------------------------------------
# Read and filter the original TSV
# ------------------------------------------------------------

message("Reading TSV:")
message("  ", input_tsv)

smorfs <- fread(
  input_tsv,
  sep = "\t",
  header = TRUE,
  na.strings = c("NA"),
  quote = "",
  data.table = TRUE,
  showProgress = TRUE
)

original_column_order <- names(smorfs)
n_input <- nrow(smorfs)

required_tsv_columns <- c(
  matching_col,
  "aa_len",
  "amp_prob",
  "flag_edge"
)

missing_tsv_columns <- setdiff(
  required_tsv_columns,
  names(smorfs)
)

if (length(missing_tsv_columns) > 0L) {
  stop(
    "The following required TSV columns are missing: ",
    paste(missing_tsv_columns, collapse = ", "),
    call. = FALSE
  )
}


# ---------------- Filter TSV ----------------

aa_len_num <- suppressWarnings(as.numeric(smorfs$aa_len))
amp_prob_num <- suppressWarnings(as.numeric(smorfs$amp_prob))
flag_edge_num <- suppressWarnings(as.numeric(smorfs$flag_edge))

keep_rows <-
  !is.na(aa_len_num) &
  aa_len_num <= max_aa_len &
  !is.na(amp_prob_num) &
  amp_prob_num > min_amp_prob &
  !is.na(flag_edge_num) &
  flag_edge_num == 0

filtered <- smorfs[keep_rows]

# Preserve the exact original order of rows after filtering.
filtered[, original_output_order__ := .I]

n_filtered <- nrow(filtered)

message("")
message("Filtering summary:")
message("  Input rows:      ", n_input)
message("  Output rows:     ", n_filtered)
message("  Rows removed:    ", n_input - n_filtered)
message("  Filters applied:")
message("    aa_len <= ", max_aa_len)
message("    amp_prob > ", min_amp_prob)
message("    flag_edge == 0")

if (n_filtered == 0L) {
  warning(
    "No rows passed the requested filters. ",
    "An output file containing only the original header will be written."
  )
  
  empty_output <- filtered[
    ,
    c(original_column_order),
    with = FALSE
  ]
  
  fwrite(
    empty_output,
    output_tsv,
    sep = "\t",
    quote = FALSE,
    na = ""
  )
  
  message("Output written:")
  message("  ", output_tsv)
  
  quit(save = "no", status = 0)
}


# Check the matching identifier.
filtered[, tsv_match_key__ := as.character(get(matching_col))]

if (anyNA(filtered$tsv_match_key__) ||
    any(!nzchar(filtered$tsv_match_key__))) {
  warning(
    "Some filtered rows have empty values in matching column '",
    matching_col,
    "'. These rows cannot match the GFF."
  )
}


# ------------------------------------------------------------
# Match against primary and secondary GFF files
# ------------------------------------------------------------

filtered_keys <- unique(
  filtered[
    !is.na(tsv_match_key__) &
      nzchar(tsv_match_key__),
    tsv_match_key__
  ]
)


# ---------------- Primary GFF ----------------

primary_gff <- read_and_parse_gff(
  input_gff,
  annotation_source = "primary_prodigal",
  match_mode = "primary"
)

primary_relevant <- primary_gff[
  gff_match_key__ %chin% filtered_keys
]

primary_duplicates <- primary_relevant[
  !is.na(gff_match_key__) &
    nzchar(gff_match_key__),
  .N,
  by = gff_match_key__
][N > 1L]

if (nrow(primary_duplicates) > 0L) {
  stop(
    "Duplicated matching keys were found in the primary GFF. ",
    "Examples: ",
    paste(
      head(primary_duplicates$gff_match_key__, 10L),
      collapse = ", "
    ),
    call. = FALSE
  )
}

primary_matched_keys <- unique(
  primary_relevant$gff_match_key__
)

unmatched_after_primary <- setdiff(
  filtered_keys,
  primary_matched_keys
)

message("")
message("Primary GFF matching:")
message("  Filtered TSV identifiers:       ", length(filtered_keys))
message("  IDs found in primary GFF:       ", length(primary_matched_keys))
message("  IDs not found in primary GFF:   ", length(unmatched_after_primary))

if (length(unmatched_after_primary) > 0L) {
  message("")
  message("  Example IDs not found in primary GFF:")
  
  for (example_id in head(unmatched_after_primary, 10L)) {
    message("    ", example_id)
  }
}


# ---------------- Secondary GFF ----------------

secondary_gff_dt <- read_and_parse_gff(
  secondary_gff,
  annotation_source = "secondary_smorfinder",
  match_mode = "exact_id"
)

secondary_relevant <- secondary_gff_dt[
  gff_match_key__ %chin% unmatched_after_primary
]

secondary_duplicates <- secondary_relevant[
  !is.na(gff_match_key__) &
    nzchar(gff_match_key__),
  .N,
  by = gff_match_key__
][N > 1L]

if (nrow(secondary_duplicates) > 0L) {
  stop(
    "Duplicated matching keys were found in the secondary GFF. ",
    "Examples: ",
    paste(
      head(secondary_duplicates$gff_match_key__, 10L),
      collapse = ", "
    ),
    call. = FALSE
  )
}

secondary_matched_keys <- unique(
  secondary_relevant$gff_match_key__
)

still_unmatched <- setdiff(
  unmatched_after_primary,
  secondary_matched_keys
)

message("")
message("Secondary GFF matching:")
message("  IDs searched in secondary GFF: ", length(unmatched_after_primary))
message("  IDs recovered from secondary:   ", length(secondary_matched_keys))
message("  IDs still unmatched:            ", length(still_unmatched))

if (length(secondary_matched_keys) > 0L) {
  message("")
  message("  Example IDs recovered from secondary GFF:")
  
  for (example_id in head(secondary_matched_keys, 10L)) {
    message("    ", example_id)
  }
}

if (length(still_unmatched) > 0L) {
  message("")
  message("  Example IDs still unmatched after both GFF files:")
  
  for (example_id in head(still_unmatched, 10L)) {
    message("    ", example_id)
  }
}


# ------------------------------------------------------------
# Combine primary and secondary annotations
# ------------------------------------------------------------

combined_gff <- rbindlist(
  list(
    primary_relevant,
    secondary_relevant
  ),
  use.names = TRUE,
  fill = TRUE
)

combined_duplicates <- combined_gff[
  !is.na(gff_match_key__) &
    nzchar(gff_match_key__),
  .N,
  by = gff_match_key__
][N > 1L]

if (nrow(combined_duplicates) > 0L) {
  stop(
    "Duplicated matching keys remained after combining GFF files.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# Prepare and append GFF columns
# ------------------------------------------------------------

gff_output_columns <- setdiff(
  names(combined_gff),
  c(
    "gff_id_final_component__",
    "gff_match_key__"
  )
)

gff_for_join <- combined_gff[
  ,
  c(
    "gff_match_key__",
    gff_output_columns
  ),
  with = FALSE
]

setkey(gff_for_join, gff_match_key__)

columns_to_append <- setdiff(
  names(gff_for_join),
  "gff_match_key__"
)

gff_for_join_copy <- copy(gff_for_join)

setnames(
  gff_for_join_copy,
  "gff_match_key__",
  "tsv_match_key__"
)

filtered <- merge(
  filtered,
  gff_for_join_copy,
  by = "tsv_match_key__",
  all.x = TRUE,
  sort = FALSE
)


# ------------------------------------------------------------
# Assess final matching success
# ------------------------------------------------------------

filtered[, gff_match_status := "unmatched"]

filtered[
  gff_annotation_source == "primary_prodigal",
  gff_match_status := "matched_primary"
]

filtered[
  gff_annotation_source == "secondary_smorfinder",
  gff_match_status := "matched_secondary"
]

n_primary <- filtered[
  gff_match_status == "matched_primary",
  .N
]

n_secondary <- filtered[
  gff_match_status == "matched_secondary",
  .N
]

n_unmatched <- filtered[
  gff_match_status == "unmatched",
  .N
]

message("")
message("Final matching summary:")
message("  Matched using primary GFF:   ", n_primary)
message("  Matched using secondary GFF: ", n_secondary)
message("  Still unmatched:             ", n_unmatched)


# ------------------------------------------------------------
# Restore output column and row order
# ------------------------------------------------------------

setorder(
  filtered,
  original_output_order__
)

# Original TSV columns always come first and remain in their
# original order.
#
# GFF fields then follow:
#   gff_seqid
#   gff_source
#   gff_type
#   gff_start
#   gff_end
#   gff_score
#   gff_strand
#   gff_phase
#   gff_attributes_raw
#   gff_attr_ID
#   gff_attr_partial
#   gff_attr_start_type
#   ...
#   gff_match_status

final_gff_columns <- c(
  gff_output_columns,
  "gff_match_status"
)

final_output_columns <- c(
  original_column_order,
  final_gff_columns
)

# Remove temporary internal columns by selecting only final columns.
result <- filtered[
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
