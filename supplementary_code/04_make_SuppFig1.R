# Internal-cohort supporting panels only; montage letters are not manuscript letters.
# A: ten melanocyte markers (Figure S1A); B: 49 RAB genes (Figure S1B).
# C: RAB7A in nine cell types (Figure S2A). These calculations are unchanged.
# Revised external S2B is generated ONLY by external_cohort/01_make_Figure_S2B.R.
# Historical external top-quartile code is deliberately not executed here.
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


# Assemble only the three internal-cohort panels.
row1 <- plot_grid(pA, pB, ncol = 2, rel_widths = c(1, 1.5),
                  labels = c("A", "B"), label_size = 14,
                  label_fontface = "bold", label_fontfamily = "Arial")
row2 <- plot_grid(pC, labels = "C", label_size = 14,
                  label_fontface = "bold", label_fontfamily = "Arial")
combined <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1, 1))
out_pdf <- file.path(fig_dir, "SuppFig_internal_scRNA_evidence.pdf")
cairo_pdf(out_pdf, width = 12, height = 10)
print(combined)
dev.off()
cat(sprintf("  Saved: %s\n", out_pdf))

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
