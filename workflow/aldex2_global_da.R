#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

log_msg <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), sprintf(...), "\n", sep = "")
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_flagstat_file <- function(path) {
  txt <- readLines(path, warn = FALSE)

  parse_metric <- function(pattern, label) {
    hit <- grep(pattern, txt, value = TRUE)
    if (length(hit) == 0) {
      return(data.table(metric = label, count = NA_real_, total = NA_real_, pct = NA_real_, line = NA_character_))
    }
    line <- hit[1]
    count <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", line)))
    pct <- suppressWarnings(as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", line)))
    total_line <- grep("in total", txt, value = TRUE)
    total <- if (length(total_line) > 0) suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", total_line[1]))) else NA_real_
    data.table(metric = label, count = count, total = total, pct = pct, line = line)
  }

  rbindlist(list(
    parse_metric(" primary mapped ", "primary_mapped"),
    parse_metric(" mapped ", "mapped")
  ), fill = TRUE)
}

collect_flagstat_summary <- function(
    stats_root = "/scratch/t.sousa/data_used/read_mapping/stats",
    envs = c("campina", "peneira", "uniforme", "riparia")) {

  all_rows <- list()
  k <- 1L

  for (env in envs) {
    env_dir <- file.path(stats_root, env)
    if (!dir.exists(env_dir)) {
      warning(sprintf("Environment directory not found, skipping: %s", env_dir))
      next
    }

    sample_dirs <- list.dirs(env_dir, recursive = FALSE, full.names = TRUE)
    for (sample_dir in sample_dirs) {
      sorted_files <- Sys.glob(file.path(sample_dir, "*.sorted.flagstat.txt"))
      if (length(sorted_files) == 0) next

      for (sorted_file in sorted_files) {
        lib_base <- basename(sorted_file)
        lib_base <- sub("\\.sorted\\.flagstat\\.txt$", "", lib_base)
        q20_file <- file.path(sample_dir, paste0(lib_base, ".primary_q20.flagstat.txt"))

        if (!file.exists(q20_file)) {
          warning(sprintf("Missing paired q20 flagstat for %s", sorted_file))
          next
        }

        dt_sorted <- read_flagstat_file(sorted_file)
        dt_sorted[, metric := fifelse(metric == "mapped", "mapped", "primary_mapped")]

        dt_q20 <- read_flagstat_file(q20_file)
        dt_q20[, metric := fifelse(metric == "mapped", "mapped_q20", "primary_mapped_q20")]

        dt <- rbind(dt_sorted, dt_q20, fill = TRUE)
        dt[, environment := env]
        dt[, sample_id := basename(sample_dir)]
        dt[, library_id := lib_base]
        dt[, sorted_flagstat := sorted_file]
        dt[, q20_flagstat := q20_file]
        all_rows[[k]] <- dt
        k <- k + 1L
      }
    }
  }

  if (length(all_rows) == 0) stop("No flagstat files found under the requested environments.", call. = FALSE)

  out <- rbindlist(all_rows, fill = TRUE)
  out[, environment := factor(environment, levels = envs)]
  out <- out[metric != "mapped_q20"]
  out[, metric := factor(metric, levels = c("mapped", "primary_mapped", "primary_mapped_q20"))]
  out[]
}

plot_flagstat_stripchart <- function(summary_dt, out_png, out_pdf = NULL, width = 9, height = 5) {
  dir_create(dirname(out_png))
  if (!is.null(out_pdf)) dir_create(dirname(out_pdf))
  
  p <- ggplot(summary_dt, aes(x = metric, y = count)) +
    geom_jitter(width = 0.18, height = 0, alpha = 0.5, size = 2) +
    
    stat_summary(
      fun = mean,
      geom = "crossbar",
      width = 0.5,
      fatten = 0,
      color = "red"
    ) +
    
    coord_cartesian(ylim = c(0, 1e7)) +
    
    facet_wrap(~ environment, nrow = 1, scales = "fixed") +
    labs(
      x = NULL,
      y = "Mapped reads (count)",
      title = "Read-mapping summary by library",
      subtitle = "Raw counts parsed from samtools flagstat outputs"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(out_png, p, width = width, height = height, dpi = 300)
  if (!is.null(out_pdf)) ggsave(out_pdf, p, width = width, height = height)
  
  invisible(p)
}

run_flagstat_step <- function(stats_root, outdir, envs) {
  dir_create(outdir)
  dt <- collect_flagstat_summary(stats_root = stats_root, envs = envs)
  fwrite(dt, file.path(outdir, "mapping_flagstat_summary.tsv"), sep = "\t")
  plot_flagstat_stripchart(
    dt,
    out_png = file.path(outdir, "mapping_flagstat_stripchart.png"),
    out_pdf = file.path(outdir, "mapping_flagstat_stripchart.pdf")
  )
  log_msg("Flagstat summary written to %s", file.path(outdir, "mapping_flagstat_summary.tsv"))
}

prepare_aldex2_input <- function(
    stats_root = "/scratch/t.sousa/data_used/read_mapping/stats",
    envs = c("campina", "peneira", "uniforme", "riparia"),
    outdir,
    remove_star = TRUE) {

  dir_create(outdir)
  idx_list <- list()
  meta_list <- list()
  k <- 1L

  for (env in envs) {
    env_dir <- file.path(stats_root, env)
    if (!dir.exists(env_dir)) {
      warning(sprintf("Environment directory not found, skipping: %s", env_dir))
      next
    }

    sample_dirs <- list.dirs(env_dir, recursive = FALSE, full.names = TRUE)
    for (sample_dir in sample_dirs) {
      idx_files <- Sys.glob(file.path(sample_dir, "*.primary_q20.idxstats.tsv"))
      if (length(idx_files) == 0) next

      for (idx in idx_files) {
        sample_name <- basename(idx)
        sample_name <- sub("\\.primary_q20\\.idxstats\\.tsv$", "", sample_name)

        dt <- fread(idx, sep = "\t", header = FALSE, col.names = c("feature_id", "feature_length", "mapped", "unmapped"))
        if (remove_star) dt <- dt[feature_id != "*"]
        dt <- dt[, .(feature_id, count = as.numeric(mapped))]
        setnames(dt, "count", sample_name)

        idx_list[[k]] <- dt
        meta_list[[k]] <- data.table(
          sample_name = sample_name,
          environment = env,
          sample_id = basename(sample_dir),
          idxstats_file = idx
        )
        k <- k + 1L
      }
    }
  }

  if (length(idx_list) == 0) stop("No *.primary_q20.idxstats.tsv files found.", call. = FALSE)

  counts_dt <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), idx_list)
  for (j in seq_along(counts_dt)) {
    if (is.numeric(counts_dt[[j]])) data.table::set(counts_dt, which(is.na(counts_dt[[j]])), j, 0)
  }

  meta_dt <- rbindlist(meta_list, fill = TRUE)
  meta_dt[, environment := factor(environment, levels = envs)]

  sample_cols <- setdiff(names(counts_dt), "feature_id")
  missing_meta <- setdiff(sample_cols, meta_dt$sample_name)
  if (length(missing_meta) > 0) stop(sprintf("Metadata missing for samples: %s", paste(missing_meta, collapse = ", ")), call. = FALSE)

  meta_dt <- meta_dt[match(sample_cols, sample_name)]

  counts_out <- file.path(outdir, "aldex2_counts_matrix.tsv")
  meta_out <- file.path(outdir, "aldex2_sample_metadata.tsv")
  fwrite(counts_dt, counts_out, sep = "\t")
  fwrite(meta_dt, meta_out, sep = "\t")

  log_msg("Counts matrix: %s", counts_out)
  log_msg("Metadata table: %s", meta_out)
  list(counts = counts_dt, metadata = meta_dt)
}

filter_low_information_cds <- function(counts_dt, min_count = 10, min_samples = 2) {
  mat <- as.matrix(counts_dt[, -1, with = FALSE])
  keep <- rowSums(mat >= min_count, na.rm = TRUE) >= min_samples
  filtered <- counts_dt[keep]
  stats <- data.table(
    n_features_before = nrow(counts_dt),
    n_features_after = nrow(filtered),
    n_removed = sum(!keep),
    min_count = min_count,
    min_samples = min_samples
  )
  list(counts = filtered, stats = stats)
}

run_aldex2_kw <- function(
    counts_tsv,
    metadata_tsv,
    outdir,
    group_col = "environment",
    mc_samples = 128,
    denom = "all",
    min_count = 10,
    min_samples = 2,
    use_mc = FALSE) {

  if (!requireNamespace("ALDEx2", quietly = TRUE)) stop("ALDEx2 is not installed in the active R environment.", call. = FALSE)

  dir_create(outdir)
  counts_dt <- fread(counts_tsv)
  meta_dt <- fread(metadata_tsv)

  if (!(group_col %in% names(meta_dt))) stop(sprintf("group_col '%s' not found in metadata.", group_col), call. = FALSE)
  if (!("sample_name" %in% names(meta_dt))) stop("Metadata must contain a 'sample_name' column.", call. = FALSE)
  if (!("feature_id" %in% names(counts_dt))) stop("Counts table must contain a 'feature_id' column.", call. = FALSE)

  sample_cols <- setdiff(names(counts_dt), "feature_id")
  if (!setequal(sample_cols, meta_dt$sample_name)) stop("Sample names in counts matrix and metadata do not match.", call. = FALSE)

  meta_dt <- meta_dt[match(sample_cols, sample_name)]
  counts_dt <- counts_dt[, c("feature_id", meta_dt$sample_name), with = FALSE]

  filt <- filter_low_information_cds(counts_dt, min_count = min_count, min_samples = min_samples)
  fwrite(filt$stats, file.path(outdir, "aldex2_filtering_summary.tsv"), sep = "\t")

  counts_filt <- filt$counts
  counts_mat <- as.matrix(counts_filt[, -1, with = FALSE])
  rownames(counts_mat) <- counts_filt$feature_id
  storage.mode(counts_mat) <- "integer"

  groups <- factor(meta_dt[[group_col]])
  names(groups) <- meta_dt$sample_name

  clr <- ALDEx2::aldex.clr(
    reads = counts_mat,
    conds = groups,
    mc.samples = mc_samples,
    denom = denom,
    useMC = use_mc,
    verbose = TRUE
  )

  kw <- ALDEx2::aldex.kw(clr, useMC = use_mc, verbose = TRUE)

  results_dt <- data.table(feature_id = rownames(kw), kw)
  raw_mean <- data.table(
    feature_id = counts_filt$feature_id,
    mean_count = rowMeans(as.matrix(counts_filt[, -1, with = FALSE])),
    prevalence = rowSums(as.matrix(counts_filt[, -1, with = FALSE]) > 0)
  )
  results_dt <- merge(results_dt, raw_mean, by = "feature_id", all.x = TRUE)
  setorder(results_dt, kw.eBH, kw.ep)

  fwrite(results_dt, file.path(outdir, "aldex2_kw_results.tsv"), sep = "\t")
  fwrite(meta_dt, file.path(outdir, "aldex2_metadata_used.tsv"), sep = "\t")
  fwrite(data.table(
    group_col = group_col,
    n_samples = ncol(counts_mat),
    n_groups = nlevels(groups),
    group_levels = paste(levels(groups), collapse = ","),
    mc_samples = mc_samples,
    denom = denom,
    min_count = min_count,
    min_samples = min_samples,
    use_mc = use_mc
  ), file.path(outdir, "aldex2_run_parameters.tsv"), sep = "\t")

  log_msg("ALDEx2 results written to %s", file.path(outdir, "aldex2_kw_results.tsv"))
  invisible(results_dt)
}

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    step = NULL,
    stats_root = "/scratch/t.sousa/data_used/read_mapping/stats",
    outdir = "results/differential_abundance/aldex2",
    envs = c("campina", "peneira", "uniforme", "riparia"),
    counts_tsv = NULL,
    metadata_tsv = NULL,
    group_col = "environment",
    mc_samples = 128,
    denom = "all",
    min_count = 10,
    min_samples = 2,
    use_mc = FALSE,
    help = FALSE
  )
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1L]] else NA_character_
    switch(key,
      "-h" = { out$help <- TRUE; i <- i + 1L },
      "--help" = { out$help <- TRUE; i <- i + 1L },
      "--step" = { out$step <- val; i <- i + 2L },
      "--stats-root" = { out$stats_root <- val; i <- i + 2L },
      "--outdir" = { out$outdir <- val; i <- i + 2L },
      "--envs" = { out$envs <- strsplit(val, ",", fixed = TRUE)[[1]]; i <- i + 2L },
      "--counts-tsv" = { out$counts_tsv <- val; i <- i + 2L },
      "--metadata-tsv" = { out$metadata_tsv <- val; i <- i + 2L },
      "--group-col" = { out$group_col <- val; i <- i + 2L },
      "--mc-samples" = { out$mc_samples <- as.integer(val); i <- i + 2L },
      "--denom" = { out$denom <- val; i <- i + 2L },
      "--min-count" = { out$min_count <- as.integer(val); i <- i + 2L },
      "--min-samples" = { out$min_samples <- as.integer(val); i <- i + 2L },
      "--use-mc" = { out$use_mc <- toupper(val) %in% c("TRUE", "T", "1"); i <- i + 2L },
      stop(sprintf("Unknown argument: %s", key), call. = FALSE)
    )
  }
  out
}

print_help <- function() {
  cat(
"Usage: Rscript workflow/aldex2_global_da.R --step <flagstat|prepare|run> [options]\n\n",
"Steps:\n",
"  --step flagstat   Parse flagstat files and generate stripchart\n",
"  --step prepare    Build ALDEx2 counts matrix + metadata from idxstats\n",
"  --step run        Run ALDEx2 Kruskal-Wallis workflow on prepared input\n"
  )
}

opt <- parse_args(args)
if (isTRUE(opt$help) || is.null(opt$step)) {
  print_help()
  quit(save = "no", status = if (isTRUE(opt$help)) 0 else 1)
}

if (opt$step == "flagstat") {
  run_flagstat_step(stats_root = opt$stats_root, outdir = opt$outdir, envs = opt$envs)
} else if (opt$step == "prepare") {
  prepare_aldex2_input(stats_root = opt$stats_root, envs = opt$envs, outdir = opt$outdir)
} else if (opt$step == "run") {
  if (is.null(opt$counts_tsv)) opt$counts_tsv <- file.path(opt$outdir, "aldex2_counts_matrix.tsv")
  if (is.null(opt$metadata_tsv)) opt$metadata_tsv <- file.path(opt$outdir, "aldex2_sample_metadata.tsv")
  run_aldex2_kw(
    counts_tsv = opt$counts_tsv,
    metadata_tsv = opt$metadata_tsv,
    outdir = opt$outdir,
    group_col = opt$group_col,
    mc_samples = opt$mc_samples,
    denom = opt$denom,
    min_count = opt$min_count,
    min_samples = opt$min_samples,
    use_mc = opt$use_mc
  )
} else {
  stop(sprintf("Unknown step: %s", opt$step), call. = FALSE)
}
