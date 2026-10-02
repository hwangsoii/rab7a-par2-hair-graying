# =============================================================================
# 06_replot_Fig1A_cluster6_umap.R — standalone keratinocyte UMAP, cluster 6
#                                    highlighted        [JID Letter Figure 1A]
#
# VISUALIZATION ONLY. This script recomputes nothing: no normalization, no
# variable-feature selection, no scaling, no PCA, no integration, no
# clustering, no UMAP, no differential expression and no KRT scoring. It reads
# keratinocytes.rds, plots the UMAP coordinates and the RNA_snn_res.0.6 cluster
# assignments that object already carries, and writes a PDF. Neither Seurat
# object is modified or re-saved.
#
# Panel spec:
#   - all 9,746 keratinocytes at their stored UMAP coordinates
#   - the 13 other clusters drawn first in light gray (#D9D9D9)
#   - cluster 6 (n = 618) drawn last in muted red (#BE3428) so it stays visible
#   - coord_fixed() to preserve the original UMAP geometry
#   - all 14 subclusters labelled by cluster number for anatomical context
#   - cluster 6 identified as the KRT35/KRT85-enriched keratinocyte cluster
#   - no legend, no axes, no ticks, no grid
#   - no panel letter, title, subtitle, cell count or explanatory sentence;
#     the letter "A" is added in Illustrator
#
# Style is inherited unchanged from config.R: Arial, theme_pub(base_size = 10).
# Page size 3.75 x 3.54 in — the same height as Fig1B_RAB7A_violin.pdf
# (05_replot_Fig1B_RAB7A_violin.R) with slightly more width for the UMAP.
#
# Input:  keratinocytes.rds, krt35_85_keratinocytes.rds
#           (both from 02_keratinocyte_analysis.R)
# Output: <figures>/Fig1A_cluster6_umap.pdf              (vector / editable)
#         <figures>/Fig1A_cluster6_umap_preview300dpi.png  (inspection only,
#                                                           NOT for submission)
# =============================================================================

source("config.R")

library(Seurat)
library(ggplot2)

require_input(c(rds_kc, rds_krt35_85),
              produced_by = "02_keratinocyte_analysis.R",
              env_hint = "HAIR_KC_DIR")

# --- Expected panel contents (asserted below, never re-derived) ---------------
EXPECTED_TOTAL     <- 9746L
EXPECTED_CLUSTERS  <- 14L
HIGHLIGHT_CLUSTER  <- "6"
EXPECTED_HIGHLIGHT <- 618L

CLUSTER_COL <- "RNA_snn_res.0.6"
REDUCTION   <- "umap"

COL_OTHER     <- "#D9D9D9"   # light gray — the 13 background clusters
COL_HIGHLIGHT <- "#BE3428"   # muted red — cluster 6 (same red as the
                             # FeaturePlots in 03_make_Fig1.R)
COL_LABEL     <- "#4D4D4D"

PT_SIZE <- 0.2               # matches DimPlot(pt.size = 0.2) in 03_make_Fig1.R

PDF_WIDTH  <- 3.75
PDF_HEIGHT <- 3.54

# =============================================================================
# Load the existing objects — read only
# =============================================================================
cat(">>> Reading", basename(rds_kc), "...\n")
keratinocytes <- readRDS(rds_kc)

# Fingerprint the object as loaded. Re-checked after plotting so that any
# accidental recomputation or mutation would be caught rather than assumed away.
state_before <- list(
  reductions = names(keratinocytes@reductions),
  assays     = names(keratinocytes@assays),
  meta_cols  = colnames(keratinocytes@meta.data),
  umap       = digest::digest(Embeddings(keratinocytes, REDUCTION)),
  clusters   = digest::digest(keratinocytes@meta.data[[CLUSTER_COL]])
)

# =============================================================================
# Verify the object before plotting anything
# =============================================================================
cat(">>> Verifying inputs ...\n")

if (ncol(keratinocytes) != EXPECTED_TOTAL) {
  stop(sprintf("keratinocytes.rds holds %d cells, expected %d. Wrong input object?",
               ncol(keratinocytes), EXPECTED_TOTAL), call. = FALSE)
}

if (!CLUSTER_COL %in% colnames(keratinocytes@meta.data)) {
  stop("keratinocytes.rds has no '", CLUSTER_COL, "' column; found: ",
       paste(colnames(keratinocytes@meta.data), collapse = ", "), call. = FALSE)
}
cluster_id <- as.character(keratinocytes@meta.data[[CLUSTER_COL]])

n_clusters <- length(unique(cluster_id))
if (n_clusters != EXPECTED_CLUSTERS) {
  stop(sprintf("%s holds %d clusters, expected %d.",
               CLUSTER_COL, n_clusters, EXPECTED_CLUSTERS), call. = FALSE)
}

is_highlight <- cluster_id == HIGHLIGHT_CLUSTER
if (sum(is_highlight) != EXPECTED_HIGHLIGHT) {
  stop(sprintf("Cluster %s holds %d cells, expected %d.",
               HIGHLIGHT_CLUSTER, sum(is_highlight), EXPECTED_HIGHLIGHT),
       call. = FALSE)
}

if (!REDUCTION %in% names(keratinocytes@reductions)) {
  stop("keratinocytes.rds has no '", REDUCTION, "' reduction; found: ",
       paste(names(keratinocytes@reductions), collapse = ", "), call. = FALSE)
}
umap <- Embeddings(keratinocytes, REDUCTION)
if (!identical(rownames(umap), colnames(keratinocytes))) {
  stop("UMAP embedding rows do not match the object's cells.", call. = FALSE)
}
if (anyNA(umap)) stop("NA values in the stored UMAP coordinates.", call. = FALSE)

# The highlighted cluster must be exactly the population carried forward into
# krt35_85_keratinocytes.rds — same cells, no more and no fewer.
cat(">>> Reading", basename(rds_krt35_85), "...\n")
krt35_85_cells <- colnames(readRDS(rds_krt35_85))

if (length(krt35_85_cells) != EXPECTED_HIGHLIGHT) {
  stop(sprintf("krt35_85_keratinocytes.rds holds %d cells, expected %d.",
               length(krt35_85_cells), EXPECTED_HIGHLIGHT), call. = FALSE)
}
highlight_cells <- colnames(keratinocytes)[is_highlight]
if (!setequal(krt35_85_cells, highlight_cells)) {
  stop(sprintf(paste("The cells in krt35_85_keratinocytes.rds are not the complete",
                     "set of cluster %s cells (%d only in the subset, %d only in",
                     "cluster %s)."),
               HIGHLIGHT_CLUSTER,
               length(setdiff(krt35_85_cells, highlight_cells)),
               length(setdiff(highlight_cells, krt35_85_cells)),
               HIGHLIGHT_CLUSTER), call. = FALSE)
}

cat(sprintf("  %d cells, %d clusters; cluster %s = %d cells == the complete krt35_85 subset.\n",
            ncol(keratinocytes), n_clusters, HIGHLIGHT_CLUSTER, sum(is_highlight)))

# =============================================================================
# Plot data — the stored coordinates, copied verbatim
#
# Split into two frames so the draw order is explicit: everything else first,
# cluster 6 last and therefore on top. Cell order within each frame is the
# object's own order, which lets the built plot be compared against the stored
# embedding element by element below.
# =============================================================================
df_all <- data.frame(
  x = umap[, 1],
  y = umap[, 2],
  cluster = cluster_id
)
df_other <- df_all[!is_highlight, c("x", "y")]
df_hl    <- df_all[ is_highlight, c("x", "y")]

# Use cluster medians as stable label anchors. This adds context without
# changing the stored UMAP or assigning unvalidated biological names.
cluster_labels <- aggregate(cbind(x, y) ~ cluster, data = df_all, FUN = median)
cluster_labels$cluster_num <- as.integer(cluster_labels$cluster)
cluster_labels <- cluster_labels[order(cluster_labels$cluster_num), ]

if (!identical(cluster_labels$cluster_num, 0:13)) {
  stop("Expected label anchors for clusters 0-13; found: ",
       paste(cluster_labels$cluster_num, collapse = ", "), call. = FALSE)
}

labels_other <- cluster_labels[cluster_labels$cluster != HIGHLIGHT_CLUSTER, ]
label_hl     <- cluster_labels[cluster_labels$cluster == HIGHLIGHT_CLUSTER, ]

# Label anchor: horizontally centred on the highlighted cluster, clear of its
# upper edge. Derived from the data, not hard-coded, and asserted to land in
# empty space so it can never sit on top of cells.
label_x   <- label_hl$x
label_y   <- max(df_hl$y) + 0.08 * diff(range(umap[, 2]))
label_pad <- c(x = 2.8, y = 1.1)   # half-extent of the callout footprint

occluded <- sum(abs(umap[, 1] - label_x) < label_pad[["x"]] &
                abs(umap[, 2] - label_y) < label_pad[["y"]])
if (occluded > 0) {
  stop(sprintf("The cluster %s callout would overlap %d cells at (%.2f, %.2f).",
               HIGHLIGHT_CLUSTER, occluded, label_x, label_y), call. = FALSE)
}

p <- ggplot() +
  geom_point(data = df_other, aes(x = x, y = y),
             colour = COL_OTHER, size = PT_SIZE, shape = 16) +
  geom_point(data = df_hl, aes(x = x, y = y),
             colour = COL_HIGHLIGHT, size = PT_SIZE, shape = 16) +
  geom_label(data = labels_other,
             aes(x = x, y = y, label = cluster),
             colour = COL_LABEL, fill = "white", label.size = NA,
             label.padding = grid::unit(0.05, "lines"),
             fontface = "bold", size = 2.35, family = "Arial") +
  geom_label(data = label_hl,
             aes(x = x, y = y, label = cluster),
             colour = COL_HIGHLIGHT, fill = "white", label.size = NA,
             label.padding = grid::unit(0.06, "lines"),
             fontface = "bold", size = 2.55, family = "Arial") +
  annotate("text", x = label_x, y = label_y,
           label = sprintf("Cluster %s\nKRT35/KRT85-enriched keratinocytes",
                           HIGHLIGHT_CLUSTER),
           colour = COL_HIGHLIGHT, fontface = "bold",
           size = 2.55, lineheight = 0.95, family = "Arial") +
  coord_fixed() +
  theme_pub() +
  theme(
    legend.position = "none",
    axis.line   = element_blank(),
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank()
  )

# =============================================================================
# Verify what the panel actually plots
#
# Coordinates and counts are read back out of the built plot, so anything
# dropped or altered between the stored embedding and the rendered layers is
# caught here rather than inferred from the input.
# =============================================================================
cat(">>> Verifying plotted cells and coordinates ...\n")
built <- ggplot_build(p)

if (length(p$layers) != 5L) {
  stop("Panel must contain two point layers and three label layers; found ",
       length(p$layers), ".", call. = FALSE)
}

n_other_plotted <- nrow(built$data[[1]])
n_hl_plotted    <- nrow(built$data[[2]])
n_total_plotted <- n_other_plotted + n_hl_plotted

cat(sprintf("  plotted: %d total (%d other clusters + %d cluster %s)\n",
            n_total_plotted, n_other_plotted, n_hl_plotted, HIGHLIGHT_CLUSTER))

if (n_total_plotted != EXPECTED_TOTAL) {
  stop(sprintf("Plotted %d cells, expected %d.", n_total_plotted, EXPECTED_TOTAL),
       call. = FALSE)
}
if (n_hl_plotted != EXPECTED_HIGHLIGHT) {
  stop(sprintf("Plotted %d cluster %s cells, expected %d.",
               n_hl_plotted, HIGHLIGHT_CLUSTER, EXPECTED_HIGHLIGHT), call. = FALSE)
}

# Every plotted coordinate must be bit-for-bit the stored coordinate.
plotted_xy <- rbind(
  cbind(built$data[[1]]$x, built$data[[1]]$y),
  cbind(built$data[[2]]$x, built$data[[2]]$y)
)
stored_xy <- rbind(umap[!is_highlight, 1:2, drop = FALSE],
                   umap[ is_highlight, 1:2, drop = FALSE])
dimnames(plotted_xy) <- NULL
dimnames(stored_xy)  <- NULL

if (!identical(plotted_xy, stored_xy)) {
  stop("Plotted UMAP coordinates differ from the stored embedding.", call. = FALSE)
}
cat("  UMAP coordinates identical to the stored embedding (bitwise).\n")

# Draw order: the highlighted cluster is the last point layer, so it renders on
# top of the gray cells.
if (!identical(unique(built$data[[1]]$colour), COL_OTHER) ||
    !identical(unique(built$data[[2]]$colour), COL_HIGHLIGHT)) {
  stop("Point colours or draw order are not as specified.", call. = FALSE)
}

# Verify cluster-number labels and the biological callout.
if (!setequal(as.character(built$data[[3]]$label), as.character(0:13)[-7])) {
  stop("Background cluster labels are incomplete or duplicated.", call. = FALSE)
}
if (!identical(as.character(built$data[[4]]$label), HIGHLIGHT_CLUSTER)) {
  stop("The highlighted cluster-number label must be '",
       HIGHLIGHT_CLUSTER, "'.", call. = FALSE)
}
expected_callout <- sprintf("Cluster %s\nKRT35/KRT85-enriched keratinocytes",
                            HIGHLIGHT_CLUSTER)
if (!identical(as.character(built$data[[5]]$label), expected_callout)) {
  stop("The biological callout is not as specified.", call. = FALSE)
}

# Nothing was recomputed and neither object was mutated.
state_after <- list(
  reductions = names(keratinocytes@reductions),
  assays     = names(keratinocytes@assays),
  meta_cols  = colnames(keratinocytes@meta.data),
  umap       = digest::digest(Embeddings(keratinocytes, REDUCTION)),
  clusters   = digest::digest(keratinocytes@meta.data[[CLUSTER_COL]])
)
if (!identical(state_before, state_after)) {
  stop("The Seurat object changed during plotting.", call. = FALSE)
}
cat(sprintf("  reductions unchanged: %s\n", paste(state_after$reductions, collapse = ", ")))
cat(sprintf("  assays unchanged:     %s\n", paste(state_after$assays, collapse = ", ")))

rm(keratinocytes); gc(verbose = FALSE)

# =============================================================================
# Export — cairo_pdf keeps text and points as vectors (editable in Illustrator)
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

out_pdf <- file.path(fig_dir, "Fig1A_cluster6_umap.pdf")
cairo_pdf(out_pdf, width = PDF_WIDTH, height = PDF_HEIGHT)
print(p)
invisible(dev.off())
cat(sprintf("  Saved: %s (%.2f x %.2f in)\n", out_pdf, PDF_WIDTH, PDF_HEIGHT))

# Raster preview for on-screen inspection only — the PDF above is the
# submission file.
out_png <- file.path(fig_dir, "Fig1A_cluster6_umap_preview300dpi.png")
png(out_png, width = PDF_WIDTH, height = PDF_HEIGHT, units = "in", res = 300)
print(p)
invisible(dev.off())
cat(sprintf("  Saved: %s (preview only, not for submission)\n", out_png))

cat("\n=== Figure 1A complete ===\n")
