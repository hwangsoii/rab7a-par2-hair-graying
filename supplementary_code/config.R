# =============================================================================
# config.R — Shared configuration for the hair graying scRNA-seq project
#
# All input and output locations are configurable. No path in this file points
# at any private or machine-specific location.
#
# Expected layout (created by the pipeline where it writes):
#
#   <HAIR_DATA_DIR>/raw/{Black1,Black2,Grey1,Grey2}/filtered_feature_bc_matrix/
#   <HAIR_DATA_DIR>/processed/                 integrated + annotated objects
#   <HAIR_DATA_DIR>/processed/keratinocyte/    keratinocyte sub-analysis objects
#   <HAIR_DATA_DIR>/public_cohort/             external-cohort objects
#   <HAIR_OUTPUT_DIR>/figures/                 figure PDFs
#   <HAIR_OUTPUT_DIR>/results/                 statistics tables
#
# Resolution order for each directory:
#   1. Its own environment variable, if set:
#        HAIR_RAW_DIR, HAIR_PROCESSED_DIR, HAIR_KC_DIR, HAIR_PUBLIC_DIR,
#        HAIR_FIGURE_DIR, HAIR_RESULTS_DIR
#   2. Derived from the roots HAIR_DATA_DIR (default <HAIR_BASE_DIR>/data) and
#      HAIR_OUTPUT_DIR (default <HAIR_BASE_DIR>), as shown above.
#
# Minimal usage:
#   HAIR_DATA_DIR=/path/to/data HAIR_OUTPUT_DIR=/path/to/output \
#     Rscript 01_analysis_pipeline.R
#
# HAIR_BASE_DIR defaults to the parent of the directory containing this file.
#
# Machine-specific paths belong in an untracked `config.local.R` beside this
# file; when present it is sourced first and may call Sys.setenv() to set any of
# the variables above. It is listed in .gitignore and is never released.
# =============================================================================

.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), mustWork = FALSE))
  }
  for (frame in rev(sys.frames())) {
    if (!is.null(frame$ofile)) {
      return(normalizePath(dirname(frame$ofile), mustWork = FALSE))
    }
  }
  normalizePath(".", mustWork = FALSE)
}

.env_dir <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else NULL
}

# Optional, untracked, machine-specific overrides (see header).
.local_config <- file.path(.get_script_dir(), "config.local.R")
if (file.exists(.local_config)) {
  source(.local_config, local = TRUE)
  cat("config.local.R sourced (untracked machine-specific paths)\n")
}

base_dir <- .env_dir("HAIR_BASE_DIR")
if (is.null(base_dir)) {
  base_dir <- normalizePath(file.path(.get_script_dir(), ".."), mustWork = FALSE)
}
if (!dir.exists(base_dir)) {
  stop("base_dir does not exist: ", base_dir,
       "\nSet HAIR_BASE_DIR to the project root, or set HAIR_DATA_DIR and ",
       "HAIR_OUTPUT_DIR directly.")
}

# --- Roots -------------------------------------------------------------------
data_dir <- .env_dir("HAIR_DATA_DIR")
if (is.null(data_dir)) data_dir <- file.path(base_dir, "data")

output_dir <- .env_dir("HAIR_OUTPUT_DIR")
if (is.null(output_dir)) output_dir <- base_dir

# --- Derived paths -----------------------------------------------------------
.resolve <- function(env_name, subdir) {
  explicit <- .env_dir(env_name)
  if (!is.null(explicit)) return(explicit)
  file.path(data_dir, subdir)
}

raw_dir        <- .resolve("HAIR_RAW_DIR",       "raw")
processed_dir  <- .resolve("HAIR_PROCESSED_DIR", "processed")
kc_dir         <- .resolve("HAIR_KC_DIR",        "processed/keratinocyte")
pub_dir        <- .resolve("HAIR_PUBLIC_DIR",    "public_cohort")

fig_dir     <- .env_dir("HAIR_FIGURE_DIR")
if (is.null(fig_dir))     fig_dir     <- file.path(output_dir, "figures")
results_dir <- .env_dir("HAIR_RESULTS_DIR")
if (is.null(results_dir)) results_dir <- file.path(output_dir, "results")

for (.d in c(processed_dir, kc_dir, fig_dir, results_dir)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

# --- RDS file paths ----------------------------------------------------------
rds_std       <- file.path(processed_dir, "data.hair.std.r0.5.rds")
rds_rename    <- file.path(processed_dir, "data.hair.rename.rds")
rds_kc        <- file.path(kc_dir, "keratinocytes.rds")
rds_krt35_85  <- file.path(kc_dir, "krt35_85_keratinocytes.rds")
rds_krt35_DEG <- file.path(kc_dir, "krt35_85_DEGs.rds")

# --- Input checks ------------------------------------------------------------
# Stop early with an actionable message instead of failing inside Read10X/readRDS.
require_input <- function(paths, produced_by = NULL, env_hint = NULL) {
  missing <- paths[!file.exists(paths) & !dir.exists(paths)]
  if (length(missing) == 0) return(invisible(TRUE))
  msg <- paste0("Missing required input(s):\n  ",
                paste(missing, collapse = "\n  "))
  if (!is.null(produced_by)) {
    msg <- paste0(msg, "\nThese are produced by: ", produced_by,
                  ". Run it first.")
  }
  if (!is.null(env_hint)) {
    msg <- paste0(msg, "\nIf they live elsewhere, set ", env_hint,
                  " to the directory that contains them.")
  }
  stop(msg, call. = FALSE)
}

require_sample_dirs <- function(samples = c("Black1", "Black2", "Grey1", "Grey2")) {
  expected <- file.path(raw_dir, samples, "filtered_feature_bc_matrix")
  require_input(expected, env_hint = "HAIR_RAW_DIR")
  invisible(expected)
}

# --- Color palettes ----------------------------------------------------------
celltype_colors <- c(
  "Endothelial"  = "#17becf",
  "Fibroblast"   = "#1f77b4",
  "Gland"        = "#ff7f0e",
  "Keratinocyte" = "#d62728",
  "Mast cell"    = "#bcbd22",
  "Melanocyte"   = "#9467bd",
  "Myeloid"      = "#8c564b",
  "Pericyte"     = "#e377c2",
  "T cell"       = "#2ca02c"
)

pub_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5"
)

hair_colors   <- c("Black" = "#1a1a2e", "Gray" = "#c8c8c8")
public_colors <- c("Black" = "#1a1a2e", "White" = "#c8c8c8")

# --- Publication labels (internal cluster names -> paper names) ---------------
pub_labels <- c(
  "Fibroblast"          = "Fibroblast",
  "Lymphocyte"          = "T cell",
  "Keratinocyte"        = "Keratinocyte",
  "Vascular"            = "Endothelial",
  "Myeloid"             = "Myeloid",
  "Pericyte"            = "Pericyte",
  "Melanocyte"          = "Melanocyte",
  "Mast cells/Basophils" = "Mast cell",
  "Gland"               = "Gland"
)

# --- Shared publication theme ------------------------------------------------
theme_pub <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") %+replace%
    ggplot2::theme(
      axis.title    = ggplot2::element_text(size = base_size),
      axis.text     = ggplot2::element_text(size = base_size - 1, color = "black"),
      legend.text   = ggplot2::element_text(size = base_size - 2),
      legend.title  = ggplot2::element_blank(),
      plot.title    = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin   = ggplot2::margin(3, 3, 3, 3)
    )
}

# --- Bonferroni helper -------------------------------------------------------
bonf_signif <- function(p_raw, n_tests) {
  p_adj <- pmin(p_raw * n_tests, 1)
  ifelse(p_adj < 0.001, "***",
         ifelse(p_adj < 0.01, "**",
                ifelse(p_adj < 0.05, "*", "ns")))
}

# --- Violin + boxplot helper -------------------------------------------------
make_violin <- function(df, y_var, y_label, colors = hair_colors) {
  fill_col <- if ("Protocol" %in% names(df)) "Protocol" else "Phenotype"
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[fill_col]], y = .data[[y_var]],
                                    fill = .data[[fill_col]])) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.4,
                         scale = "width") +
    ggplot2::geom_jitter(width = 0.15, size = 0.3, alpha = 0.15,
                         color = "black") +
    ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.9,
                          linewidth = 0.4, color = "black") +
    ggplot2::scale_fill_manual(values = colors) +
    ggpubr::stat_compare_means(method = "wilcox.test", label = "p.format",
                               label.x = 1.5, vjust = -0.5, size = 3,
                               family = "Arial") +
    ggplot2::labs(x = NULL, y = y_label) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    theme_pub(base_size = 10) +
    ggplot2::theme(legend.position = "none")
}

cat("config.R loaded.\n")
cat("  raw_dir     :", raw_dir, "\n")
cat("  processed   :", processed_dir, "\n")
cat("  kc_dir      :", kc_dir, "\n")
cat("  public_dir  :", pub_dir, "\n")
cat("  figures     :", fig_dir, "\n")
cat("  results     :", results_dir, "\n")
