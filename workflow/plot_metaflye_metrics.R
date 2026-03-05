#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript workflow/plot_metaflye_metrics.R <metrics_tsv> <plots_dir>")
}

metrics_tsv <- args[1]
plots_dir <- args[2]

if (!file.exists(metrics_tsv)) {
  stop(paste("Metrics TSV not found:", metrics_tsv))
}

dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

median_iqr <- function(x) {
  x <- x[is.finite(x)]
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  data.frame(
    y = median(x),
    ymin = qs[[1]],
    ymax = qs[[2]]
  )
}

df <- read_tsv(metrics_tsv, show_col_types = FALSE)

# Columns you want to plot (ignore non-numeric + 'site')
# We'll auto-detect numeric-ish columns safely.
numeric_cols <- names(df)[sapply(df, function(x) {
  # If already numeric/integer -> yes
  if (is.numeric(x)) return(TRUE)
  # If character, attempt numeric coercion and require at least some non-NA
  if (is.character(x)) {
    y <- suppressWarnings(as.numeric(x))
    return(sum(!is.na(y)) > 0)
  }
  FALSE
})]

# Exclude "site" explicitly if it got inferred
numeric_cols <- setdiff(numeric_cols, c("site"))

if (length(numeric_cols) == 0) {
  stop("No numeric columns found to plot.")
}

env_colors <- c(
  CAMPINA_GLOBAL  = "#FFCC00",
  UNIFORME_GLOBAL = "#99CC33",
  RIPARIA_GLOBAL  = "#3399FF",
  PENEIRA_GLOBAL  = "#FF9900"
)

for (col in numeric_cols) {
  
  # Build a per-metric plotting df with aligned value + site
  df_plot <- df %>%
    transmute(
      site = .data$site,
      val  = suppressWarnings(as.numeric(.data[[col]]))
    ) %>%
    filter(!is.na(val))
  
  if (nrow(df_plot) == 0) next
  
  p <- ggplot(df_plot, aes(x = "", y = val)) +
    
    geom_quasirandom(
      aes(fill = site),
      shape = 21,
      size = 1.5,
      dodge.width = 0.75,
      alpha = 0.8,
      color = "black"
    ) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.25,
      color = "grey40",
      alpha = 0.3
    ) +
    stat_summary(
      fun.data = median_iqr,
      geom = "errorbar",
      width = 0.15,
      color = "grey40",
      alpha = 0.3
    ) +
    
    scale_fill_manual(values = env_colors) +
    
    labs(
      title = "",
      x = NULL,
      y = col
    ) +
    
    theme_classic(base_size = 16) +
    theme(
      panel.grid       = element_blank(),
      strip.background = element_rect(color = NA),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      legend.position  = "right",
      legend.title     = element_text(size = 14),
      legend.text      = element_text(size = 12),
      axis.title       = element_text(size = 14),
      axis.text        = element_text(size = 12)
    ) +
    coord_flip()
  
  out_png <- file.path(plots_dir, paste0("stripchart_", col, ".png"))
  ggsave(out_png, plot = p, width = 6, height = 2, dpi = 200)
}

message("Plots written to: ", plots_dir)

# ggplot(
#   alpha_long,
#   aes(x = group, y = value, fill = group, color = group)
# ) +
#   geom_violin(
#     alpha = 0.25,
#     linewidth = 0,
#     position = position_dodge(width = 0.75)
#   ) +
#   geom_quasirandom(
#     shape = 21,
#     size = 1.5,
#     dodge.width = 0.75,
#     alpha = 0.5,
#     color = "black"
#   ) +
#   geom_boxplot(
#     outlier.shape = NA,
#     width = 0.15,                     
#     alpha = 0.9,
#     color = "black",
#     fill = "white"
#   ) +
#   facet_wrap(~ metric, scales = "free_y", nrow = 2) +       
#   labs(
#     x = "\nGroup",
#     y = "Alpha diversity\n",                  
#     title = ""
#   ) +
#   scale_fill_manual(name = "Group", values = env_colors) +   
#   scale_color_manual(name = "Group", values = env_colors) +  
#   theme_classic(base_size = 16) +                            
#   theme(
#     panel.grid       = element_blank(),
#     strip.background = element_rect(color = NA),
#     axis.text.x      = element_blank(),     
#     axis.ticks.x     = element_blank(),     
#     legend.position  = "right",
#     legend.title     = element_text(size = 14),
#     legend.text      = element_text(size = 12),
#     axis.title       = element_text(size = 14),
#     axis.text        = element_text(size = 12)
#   )
