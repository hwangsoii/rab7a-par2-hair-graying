# =============================================================================
# 04_make_SuppFig1.R — Supplementary scRNA-seq supporting evidence montage
#
# As in 03_make_Fig1.R, the montage letters are internal identifiers. The
# published supplementary figures are assembled in Illustrator:
#
#   montage pA  melanocyte gene lollipop (Bonferroni)   -> Figure S1A
#   montage pB  RAB GTPase family screen (Bonferroni)   -> Figure S1B
#   montage pC  RAB7A across all cell types             -> Figure S2A
#   montage pD  RAB7A in KRT35/85+ KC, external cohort  -> Figure S2B
#
# Manuscript Figure S1C-S1E (KRT31/KRT35/KRT85 feature plots) are produced by
# 03_make_Fig1.R.
#
# Panel D reads the exact normalized source-data values provided in
# ../data/external/Figure_S2B_source_data.tsv. The table contains cells from
# the top quartile of public-cohort keratinocyte subclusters ranked by mean
# KRT35/KRT85 expression; see supplementary_code/external_cohort/README.md.
#
# Input:  data.hair.std.r0.5.rds, data.hair.rename.rds,
#         krt35_85_keratinocytes.rds, Figure_S2B_source_data.tsv
# Output: <figures>/SuppFig_scRNA_evidence_v2.pdf,
#         <results>/SuppFig_scRNA_evidence_stats.csv
# =============================================================================

source("config.R")

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(cowplot)

require_input(c(rds_std, rds_rename, rds_krt35_85),
              produced_by = "01_analysis_pipeline.R and 02_keratinocyte_analysis.R",
              env_hint = "HAIR_PROCESSED_DIR / HAIR_KC_DIR")

stats_list <- list()

# =============================================================================
# PANEL A — Melanocyte markers NOT different (lollipop)  [Figure S1A]
# =============================================================================
cat(">>> Panel A: Melanocyte markers ...\n")

hair.std <- readRDS(rds_std)
MC <- subset(hair.std, idents = c("13"))
DefaultAssay(MC) <- "RNA"
MC <- JoinLayers(MC)

# Key melanocyte genes: melanogenesis enzymes, identity, transport, signaling
mc_genes <- c("TYR", "TYRP1", "DCT", "MLANA", "PMEL", "MITF",
              "KIT", "MC1R", "RAB27A", "MYO5A")
expr_mat <- GetAssayData(MC, slot = "data")
mc_prot  <- MC$protocol

lollipop_rows <- list()
for (g in mc_genes) {
  if (g %in% rownames(expr_mat)) {
    vals <- expr_mat[g, ]
    b_vals <- vals[mc_prot == "Black"]
    g_vals <- vals[mc_prot == "Grey"]
    wt <- wilcox.test(b_vals, g_vals)
    log2fc <- log2(mean(expm1(b_vals)) + 1) - log2(mean(expm1(g_vals)) + 1)
    lollipop_rows[[length(lollipop_rows) + 1]] <- data.frame(
      Gene = g, log2FC = log2fc, wilcox_p = wt$p.value)
  }
}
rm(MC); gc()

df_lollipop <- bind_rows(lollipop_rows)
n_mc_genes <- nrow(df_lollipop)
df_lollipop$p_bonf <- pmin(df_lollipop$wilcox_p * n_mc_genes, 1)
df_lollipop$sig_label <- bonf_signif(df_lollipop$wilcox_p, n_mc_genes)
df_lollipop$Gene <- factor(df_lollipop$Gene, levels = rev(mc_genes))

pA <- ggplot(df_lollipop, aes(x = log2FC, y = Gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2FC, y = Gene, yend = Gene),
               linewidth = 0.8, color = "grey60") +
  geom_point(size = 3.5, color = "grey60") +
  geom_text(aes(label = sig_label), hjust = -0.5, vjust = 0.3, size = 3.5,
            family = "Arial", color = "grey40") +
  labs(x = expression(Log[2]~"FC (Black vs Gray, Wilcoxon, Bonferroni)"),
       y = NULL) +
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  theme_pub(base_size = 10) +
  theme(axis.text.y = element_text(face = "italic"))

# =============================================================================
# PANEL B — RAB GTPase family unbiased screen (Bonferroni)  [Figure S1B]
# =============================================================================
cat(">>> Panel B: RAB GTPase family ...\n")
krt35_85 <- readRDS(rds_krt35_85)
DefaultAssay(krt35_85) <- "RNA"
krt35_85 <- JoinLayers(krt35_85)

all_genes <- rownames(krt35_85)
all_rab_genes <- sort(all_genes[grepl("^RAB[0-9]", all_genes)])

prot <- factor(krt35_85$protocol, levels = c("Black", "Grey"))
rab_results <- data.frame()

for (g in all_rab_genes) {
  expr <- GetAssayData(krt35_85, slot = "data")[g, ]
  pct <- mean(expr > 0) * 100
  if (pct > 5) {
    b_vals <- expr[prot == "Black"]
    g_vals <- expr[prot == "Grey"]
    wt <- wilcox.test(b_vals, g_vals)
    log2fc <- log2(mean(expm1(b_vals)) + 1) - log2(mean(expm1(g_vals)) + 1)
    rab_results <- rbind(rab_results, data.frame(
      gene = g, pct_expressing = pct,
      mean_b = mean(b_vals), mean_g = mean(g_vals),
      log2fc = log2fc, wilcox_p = wt$p.value, stringsAsFactors = FALSE))
  }
}
rab_results <- rab_results %>% arrange(wilcox_p)
rab_results$rank <- seq_len(nrow(rab_results))
n_total_tested <- nrow(rab_results)
rab_results$p_bonf <- pmin(rab_results$wilcox_p * n_total_tested, 1)
rab_results$signif_bonf <- bonf_signif(rab_results$wilcox_p, n_total_tested)

cat(sprintf("  %d RAB genes tested, RAB7A rank: %d\n",
            n_total_tested, rab_results$rank[rab_results$gene == "RAB7A"]))

top_n <- 10
top_rab <- rab_results[1:top_n, ]
top_rab$gene_f <- factor(top_rab$gene, levels = rev(top_rab$gene))
rab_colors <- ifelse(top_rab$gene == "RAB7A", "#d62728",
              ifelse(top_rab$p_bonf < 0.05, "#1f77b4", "grey60"))
names(rab_colors) <- top_rab$gene
gene_faces_b <- ifelse(rev(top_rab$gene) == "RAB7A", "bold.italic", "italic")

pB <- ggplot(top_rab, aes(x = log2fc, y = gene_f)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted",
             color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2fc, y = gene_f, yend = gene_f, color = gene),
               linewidth = 0.6, show.legend = FALSE) +
  geom_point(aes(color = gene), size = 2.5, show.legend = FALSE) +
  geom_text(aes(label = signif_bonf), hjust = -0.3, vjust = 0.3, size = 2.5,
            family = "Arial", color = "grey40") +
  scale_color_manual(values = rab_colors) +
  labs(x = expression(Log[2]~"FC (Black vs Gray, Wilcoxon, Bonferroni)"),
       y = NULL,
       subtitle = sprintf("Top 10 of %d RAB GTPases in KRT35/85+ KC", n_total_tested)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_pub(base_size = 9) +
  theme(axis.text.y = element_text(face = gene_faces_b, size = 7),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40",
                                     family = "Arial"))

# Stats for Panel B
prot_krt <- krt35_85$protocol
for (i in seq_len(nrow(rab_results))) {
  stats_list[[length(stats_list) + 1]] <- data.frame(
    panel = "B", cell_type = paste0("KRT35/85+ KC (", rab_results$gene[i], ")"),
    cohort = "internal",
    n_black = sum(prot_krt == "Black"), n_grey = sum(prot_krt == "Grey"),
    mean_black = rab_results$mean_b[i], mean_grey = rab_results$mean_g[i],
    fold_change = rab_results$mean_b[i] / max(rab_results$mean_g[i], 1e-10),
    wilcox_p = rab_results$wilcox_p[i],
    bonferroni_p = rab_results$p_bonf[i],
    significant = rab_results$p_bonf[i] < 0.05)
}
rm(krt35_85); gc()

# =============================================================================
# PANEL C — RAB7A across ALL cell types (internal, Bonferroni)  [Figure S2A]
# =============================================================================
cat(">>> Panel C: RAB7A across all cell types ...\n")
rename <- readRDS(rds_rename)
DefaultAssay(rename) <- "RNA"

rename$cell_type_pub <- unname(pub_labels[as.character(Idents(rename))])
rename$Protocol <- factor(rename$protocol, levels = c("Black", "Grey"),
                          labels = c("Black", "Gray"))

df_all <- data.frame(
  RAB7A    = FetchData(rename, vars = "RAB7A")[, 1],
  CellType = rename$cell_type_pub,
  Protocol = rename$Protocol)

ct_order <- df_all %>%
  group_by(CellType) %>%
  summarise(diff = mean(RAB7A[Protocol == "Black"]) -
                   mean(RAB7A[Protocol == "Gray"]), .groups = "drop") %>%
  arrange(desc(diff)) %>% pull(CellType)
df_all$CellType <- factor(df_all$CellType, levels = ct_order)

n_celltypes <- length(ct_order)
ct_stats <- df_all %>%
  group_by(CellType) %>%
  summarise(p_raw = wilcox.test(RAB7A[Protocol == "Black"],
                                RAB7A[Protocol == "Gray"])$p.value, .groups = "drop")
ct_stats$signif_bonf <- bonf_signif(ct_stats$p_raw, n_celltypes)

ct_counts <- df_all %>% group_by(CellType) %>% summarise(n = n(), .groups = "drop")
ct_labels <- setNames(
  paste0(ct_counts$CellType, "\n(n=", format(ct_counts$n, big.mark = ","), ")"),
  ct_counts$CellType)

y_max_ct <- max(df_all$RAB7A, na.rm = TRUE)
annot_ct <- data.frame(
  CellType = factor(ct_stats$CellType, levels = ct_order),
  label = ct_stats$signif_bonf, y = y_max_ct * 1.08)

pC <- ggplot(df_all, aes(x = CellType, y = RAB7A, fill = Protocol)) +
  geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.3,
              position = position_dodge(0.8), scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.9, linewidth = 0.3,
               position = position_dodge(0.8), color = "black") +
  scale_fill_manual(values = hair_colors) +
  geom_text(data = annot_ct, aes(x = CellType, y = y, label = label),
            inherit.aes = FALSE, size = 3, family = "Arial") +
  scale_x_discrete(labels = ct_labels) +
  labs(x = NULL, y = "RAB7A expression",
       subtitle = sprintf("Bonferroni corrected (%d cell types)", n_celltypes)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "top", legend.key.size = unit(0.3, "cm"),
        plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey40",
                                     family = "Arial"))

for (ct in ct_order) {
  sub <- df_all %>% filter(CellType == ct)
  b <- sub %>% filter(Protocol == "Black") %>% pull(RAB7A)
  g <- sub %>% filter(Protocol == "Gray") %>% pull(RAB7A)
  wt <- wilcox.test(b, g)
  stats_list[[length(stats_list) + 1]] <- data.frame(
    panel = "C", cell_type = ct, cohort = "internal",
    n_black = length(b), n_grey = length(g),
    mean_black = mean(b), mean_grey = mean(g),
    fold_change = mean(b) / max(mean(g), 1e-10),
    wilcox_p = wt$p.value,
    bonferroni_p = min(wt$p.value * n_celltypes, 1),
    significant = min(wt$p.value * n_celltypes, 1) < 0.05)
}
rm(hair.std, rename); gc()

# =============================================================================
# PANEL D — RAB7A in KRT35/85-enriched KC (external cohort)  [Figure S2B]
# =============================================================================
cat(">>> Panel D: RAB7A in KRT35/85+ KC (public) ...\n")

external_s2b_path <- Sys.getenv(
  "HAIR_EXTERNAL_S2B_SOURCE",
  unset = file.path("..", "data", "external", "Figure_S2B_source_data.tsv")
)
if (file.exists(external_s2b_path)) {
  external_s2b <- read.delim(external_s2b_path, stringsAsFactors = FALSE)
  required_external_columns <- c("phenotype", "RAB7A")
  missing_external_columns <- setdiff(required_external_columns, names(external_s2b))
  if (length(missing_external_columns) > 0) {
    stop("Figure S2B source data are missing columns: ",
         paste(missing_external_columns, collapse = ", "))
  }

  df_d <- data.frame(
    RAB7A = external_s2b$RAB7A,
    Phenotype = external_s2b$phenotype
  ) %>%
    filter(!is.na(Phenotype), !is.na(RAB7A))
  df_d$Phenotype <- factor(df_d$Phenotype, levels = c("black", "white"),
                           labels = c("Black", "White"))

  pD <- make_violin(df_d, "RAB7A", "RAB7A expression", colors = public_colors) +
    ggtitle("KRT35/85+ KC (public)") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                    family = "Arial"))

  b <- df_d %>% filter(Phenotype == "Black") %>% pull(RAB7A)
  w <- df_d %>% filter(Phenotype == "White") %>% pull(RAB7A)
  wt <- wilcox.test(b, w)
  stats_list[[length(stats_list) + 1]] <- data.frame(
    panel = "D", cell_type = "KRT35/85+ KC", cohort = "public_v2",
    n_black = length(b), n_grey = length(w),
    mean_black = mean(b), mean_grey = mean(w),
    fold_change = mean(b) / max(mean(w), 1e-10),
    wilcox_p = wt$p.value, bonferroni_p = NA, significant = wt$p.value < 0.05)

  rm(external_s2b); gc()
} else {
  cat("  WARNING: Figure S2B source data not found, skipping Panel D\n")
  pD <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "Panel D\n(public data not available)")
}

# =============================================================================
# ASSEMBLE 4-PANEL FIGURE
# =============================================================================
cat(">>> Assembling supplementary figure ...\n")

row1 <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1.5),
                  labels = c("A", "B"), label_size = 14,
                  label_fontface = "bold", label_fontfamily = "Arial")

row2 <- plot_grid(pC, pD, ncol = 2, rel_widths = c(2, 1),
                  labels = c("C", "D"), label_size = 14,
                  label_fontface = "bold", label_fontfamily = "Arial")

combined <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))

out_pdf <- file.path(fig_dir, "SuppFig_scRNA_evidence_v2.pdf")
cairo_pdf(out_pdf, width = 12, height = 10)
print(combined)
dev.off()
cat(sprintf("  Saved: %s\n", out_pdf))

# =============================================================================
# SAVE STATISTICS
# =============================================================================
stats_df <- bind_rows(stats_list)
stats_csv <- file.path(results_dir, "SuppFig_scRNA_evidence_stats.csv")
write.csv(stats_df, stats_csv, row.names = FALSE)
cat(sprintf("  Saved: %s\n", stats_csv))

cat("\n=== Summary ===\n")
print(stats_df %>% dplyr::select(panel, cell_type, cohort,
                                  wilcox_p, bonferroni_p, significant))
cat("\n=== Supplementary Figure complete ===\n")
