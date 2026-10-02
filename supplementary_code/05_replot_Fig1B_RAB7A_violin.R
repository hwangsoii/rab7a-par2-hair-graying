# =============================================================================
# 05_replot_Fig1B_RAB7A_violin.R — standalone replot of the RAB7A violin panel
#
# Manuscript Figure 1B (this panel is generated as `pH` inside the
# 03_make_Fig1.R montage; the montage letters are internal identifiers only).
#
# The default run REPLOTS ONLY. --nominal-cell-p also computes the explicitly
# requested pooled-cell Wilcoxon test in a separate output. No clustering or
# normalization is performed: the script reads the already-saved
# krt35_85_keratinocytes.rds and plots the normalized RAB7A values that object
# already carries.
#
# Panel spec:
#   - all 618 cells of the previously selected keratinocyte cluster 6
#     (RNA_snn_res.0.6 == 6), Black n = 382 / Gray n = 236
#   - violin with embedded boxplot, normalized RAB7A expression
#   - default: no P annotation; --nominal-cell-p: nominal P, methods in legend
#   - no donor-level dots, donor means or paired lines
#   - two-line title linking the comparison directly to cluster 6 in panel A
#   - panel letter and legend are added in Illustrator
#
# Style is inherited unchanged from config.R: Arial, theme_pub(base_size = 10),
# hair_colors (Black #1a1a2e / Gray #c8c8c8). The page size is the footprint
# this panel occupies in Fig1_combined_BtoH.pdf, so the output drops into the
# Illustrator layout at 100% scale.
#
# Input:  krt35_85_keratinocytes.rds   (from 02_keratinocyte_analysis.R)
# Output: <figures>/Fig1B_RAB7A_violin.pdf   (vector / editable)
# =============================================================================

source("config.R")

library(Seurat)
library(ggplot2)
nominal_cell_p <- "--nominal-cell-p" %in% commandArgs(trailingOnly = TRUE)
input_sha <- digest::digest(file = rds_krt35_85, algo = "sha256", serialize = FALSE)

require_input(rds_krt35_85,
              produced_by = "02_keratinocyte_analysis.R",
              env_hint = "HAIR_KC_DIR")

# --- Expected panel contents (asserted below, never re-derived) ---------------
EXPECTED_TOTAL <- 618L
EXPECTED_N     <- c(Black = 382L, Gray = 236L)

# Footprint of this panel inside Fig1_combined_BtoH.pdf:
#   width  = 8.27 * 3.5 / (6.5 + 0.5 + 3.5)
#   height = 9.00 * 5.5 / (4.5 + 4.0 + 5.5)
PDF_WIDTH  <- 2.76
PDF_HEIGHT <- 3.54

# =============================================================================
# Load the existing object — read only
# =============================================================================
cat(">>> Reading", basename(rds_krt35_85), "...\n")
krt35_85 <- readRDS(rds_krt35_85)
DefaultAssay(krt35_85) <- "RNA"

if (ncol(krt35_85) != EXPECTED_TOTAL) {
  stop(sprintf("Object holds %d cells, expected %d. Wrong input object?",
               ncol(krt35_85), EXPECTED_TOTAL), call. = FALSE)
}

# All cells come from keratinocyte cluster 6; confirm before plotting "all 618".
parent_clusters <- unique(as.character(krt35_85$RNA_snn_res.0.6))
if (!identical(parent_clusters, "6")) {
  stop("Expected every cell to come from keratinocyte cluster 6, found: ",
       paste(parent_clusters, collapse = ", "), call. = FALSE)
}
cat(sprintf("  %d cells, all from keratinocyte cluster 6.\n", ncol(krt35_85)))

# =============================================================================
# Pull normalized RAB7A straight from the RNA data layer
#
# The object also carries a RAB7A meta.data column written by
# 02_keratinocyte_analysis.R; it is asserted to be identical to the data layer
# so the panel cannot silently plot a stale copy.
# =============================================================================
rab7a <- as.numeric(LayerData(krt35_85[["RNA"]], layer = "data")["RAB7A", ])

if (!is.null(krt35_85$RAB7A) &&
    !isTRUE(all.equal(unname(krt35_85$RAB7A), rab7a))) {
  stop("RAB7A meta.data column disagrees with the RNA data layer.", call. = FALSE)
}
if (anyNA(rab7a)) stop("NA values in normalized RAB7A.", call. = FALSE)

rab7a_data <- data.frame(
  RAB7A    = rab7a,
  Protocol = factor(krt35_85$protocol,
                    levels = c("Black", "Grey"),
                    labels = c("Black", "Gray"))
)

if (anyNA(rab7a_data$Protocol)) {
  stop("Unmapped protocol label(s): ",
       paste(setdiff(unique(krt35_85$protocol), c("Black", "Grey")),
             collapse = ", "), call. = FALSE)
}

cell_test <- NULL
if (nominal_cell_p) {
  black <- rab7a_data$RAB7A[rab7a_data$Protocol == "Black"]
  gray <- rab7a_data$RAB7A[rab7a_data$Protocol == "Gray"]
  result <- stats::wilcox.test(black, gray, alternative = "two.sided",
                              paired = FALSE, exact = FALSE, correct = TRUE)
  cell_test <- list(
    method = "Wilcoxon rank-sum", alternative = "two.sided", paired = FALSE,
    exact = FALSE, continuity_correction = TRUE, ties_adjusted = TRUE,
    n_black = length(black), n_gray = length(gray),
    zeros_black = sum(black == 0), zeros_gray = sum(gray == 0),
    statistic_w = unname(result$statistic), p_nominal = result$p.value,
    p_label = paste0("Nominal P = ", format.pval(result$p.value, digits = 3, eps = 1e-300)),
    zeros_retained = TRUE, donor_adjusted = FALSE, multiple_testing_adjusted = FALSE
  )
  stopifnot(is.finite(result$p.value), result$p.value > 0, result$p.value <= 1)
}

rm(krt35_85); gc(verbose = FALSE)

# =============================================================================
# Panel
# =============================================================================
p <- ggplot(rab7a_data, aes(x = Protocol, y = RAB7A, fill = Protocol)) +
  geom_violin(trim = FALSE, alpha = 0.7, linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9, linewidth = 0.4) +
  scale_fill_manual(values = hair_colors) +
  labs(
    title = "Cluster 6\nKRT35/KRT85-enriched keratinocytes",
    x = NULL,
    y = "Normalized RAB7A expression"
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      family = "Arial", face = "bold", size = 8.5,
      hjust = 0.5, lineheight = 0.95,
      margin = margin(b = 3)
    )
  )

# =============================================================================
# Verify what the panel actually plots
#
# Counts are taken from the built plot, not from the input frame, so anything
# dropped between the data and the rendered layers would be caught here.
# =============================================================================
cat(">>> Verifying plotted cell totals ...\n")
built <- ggplot_build(p)

n_input <- table(rab7a_data$Protocol)

violin_layer <- built$data[[1]]
n_violin <- if ("n" %in% names(violin_layer)) {
  stats::setNames(
    vapply(split(violin_layer$n, violin_layer$group), function(v) as.integer(unique(v)[1]), integer(1)),
    levels(rab7a_data$Protocol)
  )
} else {
  NULL   # older ggplot2 does not expose per-group n from stat_ydensity
}

report <- function(label, counts) {
  cat(sprintf("  %-22s Black = %d, Gray = %d, total = %d\n",
              label, counts[["Black"]], counts[["Gray"]], sum(counts)))
}
report("input data:", n_input)
if (!is.null(n_violin)) report("violin layer (built):", n_violin)

check <- function(counts, what) {
  if (!identical(as.integer(counts[c("Black", "Gray")]), as.integer(EXPECTED_N))) {
    stop(sprintf("%s cell totals are Black = %d / Gray = %d, expected %d / %d.",
                 what, counts[["Black"]], counts[["Gray"]],
                 EXPECTED_N[["Black"]], EXPECTED_N[["Gray"]]), call. = FALSE)
  }
}
check(n_input, "Input")
if (!is.null(n_violin)) check(n_violin, "Plotted (violin layer)")
if (sum(n_input) != EXPECTED_TOTAL) {
  stop("Plotted total is not ", EXPECTED_TOTAL, ".", call. = FALSE)
}
cat("  OK: 382 Black / 236 Gray / 618 total.\n")

# The boxplot layer must summarise the same cells (no outlier trimming of data).
box_layer <- built$data[[2]]
if (nrow(box_layer) != 2L) {
  stop("Boxplot layer has ", nrow(box_layer), " groups, expected 2.", call. = FALSE)
}

# Guard the base geometry before adding the requested optional annotation.
layer_stats <- vapply(p$layers, function(l) class(l$stat)[1], character(1))
if (length(p$layers) != 2L ||
    !identical(unname(layer_stats), c("StatYdensity", "StatBoxplot"))) {
  stop("Panel must contain exactly one violin and one boxplot layer; found: ",
       paste(layer_stats, collapse = ", "), call. = FALSE)
}

if (nominal_cell_p) {
  p <- p + annotate("text", x = 1.5, y = Inf, vjust = 1,
                    label = cell_test$p_label, family = "Arial", size = 2.8) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  annotated <- ggplot_build(p)
  stopifnot(identical(built$data[1:2], annotated$data[1:2]))
}

# =============================================================================
# Export — cairo_pdf keeps text and shapes as vectors (editable in Illustrator)
# =============================================================================
if (requireNamespace("systemfonts", quietly = TRUE)) {
  # match_fonts() in systemfonts >= 1.1.0, match_font() before that.
  arial_path <- if (exists("match_fonts", asNamespace("systemfonts"))) {
    systemfonts::match_fonts("Arial")$path
  } else {
    systemfonts::match_font("Arial")$path
  }
  if (!grepl("arial", arial_path, ignore.case = TRUE)) {
    warning("Arial did not resolve to an Arial font file (got: ", arial_path,
            "); check the output PDF.", call. = FALSE)
  }
}

out_pdf <- file.path(fig_dir, if (nominal_cell_p) "Fig1B_RAB7A_violin_nominalP.pdf"
                    else "Fig1B_RAB7A_violin.pdf")
if (nominal_cell_p && file.exists(out_pdf)) {
  stop("Nominal-P output already exists; preserve it before regenerating.", call. = FALSE)
}
cairo_pdf(out_pdf, width = PDF_WIDTH, height = PDF_HEIGHT)
print(p)
invisible(dev.off())

stopifnot(identical(input_sha, digest::digest(file = rds_krt35_85, algo = "sha256", serialize = FALSE)))
if (nominal_cell_p) {
  preview <- sub("\\.pdf$", ".png", out_pdf)
  png(preview, width = PDF_WIDTH, height = PDF_HEIGHT, units = "in", res = 300,
      type = "cairo", bg = "white")
  print(p)
  invisible(dev.off())
  relative <- function(path) {
    prefix <- paste0(normalizePath(base_dir), "/")
    full <- normalizePath(path)
    if (startsWith(full, prefix)) return(substring(full, nchar(prefix) + 1L))
    paste0("external/", basename(full))
  }
  file_record <- function(path) list(path = relative(path),
    sha256 = digest::digest(file = path, algo = "sha256", serialize = FALSE))
  qa <- list(
    source = file_record(rds_krt35_85),
    script = file_record(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
    cluster = "6", total_cells = EXPECTED_TOTAL,
    plotted_n = as.list(n_violin), cell_test = cell_test,
    base_geometry_unchanged = TRUE,
    test_method_location = "Figure legend and manuscript Methods, not inside the panel",
    design = "Four libraries from two paired donors; cells are not independent donor replicates",
    outputs = lapply(c(out_pdf, preview), file_record)
  )
  jsonlite::write_json(qa, file.path(fig_dir, "Fig1B_RAB7A_violin_nominalP_QA.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat("Nominal cell-level P:", format(result$p.value, digits = 16), "\n")
}

cat(sprintf("  Saved: %s (%.2f x %.2f in)\n", out_pdf, PDF_WIDTH, PDF_HEIGHT))
cat("\n=== Figure 1B replot complete ===\n")
