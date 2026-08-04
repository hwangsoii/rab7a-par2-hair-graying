# =============================================================================
# 02_keratinocyte_analysis.R — KC subclustering, KRT35/85 identification, DEGs
#
# Consolidates fig2.r + KRT35_85_DEG_Analysis.r into one clean script.
#
# Input:  data.hair.std.r0.5.rds (from 01_analysis_pipeline.R)
# Output: keratinocytes.rds, krt35_85_keratinocytes.rds, krt35_85_DEGs.rds
# =============================================================================

source("config.R")

library(Seurat)
library(dplyr)
library(ggplot2)

# =============================================================================
# STEP 1: Subset keratinocytes and re-cluster
# =============================================================================
cat(">>> Loading integrated object ...\n")
require_input(rds_std, produced_by = "01_analysis_pipeline.R",
              env_hint = "HAIR_PROCESSED_DIR")
hair.std.r0.5 <- readRDS(rds_std)

keratinocyte_clusters <- c(2, 3, 5, 8, 10, 11, 12)
keratinocytes <- subset(hair.std.r0.5, subset = seurat_clusters %in% keratinocyte_clusters)
rm(hair.std.r0.5); gc()

DefaultAssay(keratinocytes) <- "RNA"
keratinocytes <- NormalizeData(keratinocytes)
keratinocytes <- FindVariableFeatures(keratinocytes, selection.method = "vst", nfeatures = 2000)
keratinocytes <- ScaleData(keratinocytes)
keratinocytes <- RunPCA(keratinocytes, npcs = 30)
keratinocytes <- FindNeighbors(keratinocytes, reduction = "pca", dims = 1:20)
keratinocytes <- FindClusters(keratinocytes, resolution = 0.6)
keratinocytes <- RunUMAP(keratinocytes, reduction = "pca", dims = 1:20,
                         n.neighbors = 30, min.dist = 0.1, spread = 1.2)
keratinocytes <- JoinLayers(keratinocytes)

cat(sprintf("  Keratinocyte cells: %d, clusters: %d\n",
            ncol(keratinocytes), length(unique(keratinocytes$seurat_clusters))))

# =============================================================================
# STEP 2: Identify KRT35/85-enriched clusters
# =============================================================================
cat(">>> Identifying KRT35/85 clusters ...\n")

for (gene in c("KRT35", "KRT85")) {
  keratinocytes[[gene]] <- FetchData(keratinocytes, vars = gene)[[gene]]
}
keratinocytes$KRT35_KRT85_score <- (keratinocytes$KRT35 + keratinocytes$KRT85) / 2

cluster_krt35_85 <- keratinocytes@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_krt35_85 = mean(KRT35_KRT85_score),
    total_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_krt35_85))

krt35_85_threshold <- 0.5
krt35_85_clusters <- cluster_krt35_85$seurat_clusters[cluster_krt35_85$mean_krt35_85 > krt35_85_threshold]

cat("  KRT35/85 expression by cluster:\n")
print(as.data.frame(cluster_krt35_85))
cat(sprintf("  High KRT35/85 clusters (mean > %.1f): %s\n",
            krt35_85_threshold, paste(krt35_85_clusters, collapse = ", ")))

# =============================================================================
# STEP 3: Subset KRT35/85 clusters, re-cluster, run DEG
# =============================================================================
cat(">>> Subsetting KRT35/85 clusters and re-clustering ...\n")

krt35_85_kc <- subset(keratinocytes, subset = seurat_clusters %in% krt35_85_clusters)
DefaultAssay(krt35_85_kc) <- "RNA"
krt35_85_kc <- NormalizeData(krt35_85_kc)
krt35_85_kc <- FindVariableFeatures(krt35_85_kc, nfeatures = 1500)
krt35_85_kc <- ScaleData(krt35_85_kc)
krt35_85_kc <- RunPCA(krt35_85_kc, npcs = 20)
krt35_85_kc <- FindNeighbors(krt35_85_kc, reduction = "pca", dims = 1:15)
krt35_85_kc <- FindClusters(krt35_85_kc, resolution = 0.5)
krt35_85_kc <- RunUMAP(krt35_85_kc, reduction = "pca", dims = 1:15)

for (gene in c("KRT31", "KRT35", "KRT85", "RAB7A")) {
  krt35_85_kc[[gene]] <- FetchData(krt35_85_kc, vars = gene)[[gene]]
}

cat(sprintf("  KRT35/85 KC cells: %d (Black: %d, Grey: %d)\n",
            ncol(krt35_85_kc),
            sum(krt35_85_kc$protocol == "Black"),
            sum(krt35_85_kc$protocol == "Grey")))

# DEG: Black vs Grey in KRT35/85+ keratinocytes
cat(">>> Running DEG analysis (Black vs Grey) ...\n")
Idents(krt35_85_kc) <- "protocol"
krt35_85_DEGs <- FindMarkers(krt35_85_kc,
                              ident.1 = "Black", ident.2 = "Grey",
                              min.pct = 0.1, logfc.threshold = 0.25,
                              test.use = "wilcox")
krt35_85_DEGs$gene <- rownames(krt35_85_DEGs)

n_sig <- sum(krt35_85_DEGs$p_val_adj < 0.05)
cat(sprintf("  DEGs found: %d total, %d significant (adj p < 0.05)\n",
            nrow(krt35_85_DEGs), n_sig))

# Check RAB7A rank
if ("RAB7A" %in% krt35_85_DEGs$gene) {
  rab7a_row <- krt35_85_DEGs["RAB7A", ]
  cat(sprintf("  RAB7A: log2FC = %.3f, p_val_adj = %.2e\n",
              rab7a_row$avg_log2FC, rab7a_row$p_val_adj))
}

# =============================================================================
# STEP 4: Identify KRT31 clusters
# =============================================================================
cluster_krt31 <- krt35_85_kc@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_krt31 = mean(KRT31), .groups = "drop") %>%
  arrange(desc(mean_krt31))

krt31_threshold <- 0.5
krt31_clusters <- cluster_krt31$seurat_clusters[cluster_krt31$mean_krt31 > krt31_threshold]
cat(sprintf("  KRT31-high clusters: %s\n", paste(krt31_clusters, collapse = ", ")))

# =============================================================================
# STEP 5: Save outputs
# =============================================================================
cat(">>> Saving RDS files ...\n")

saveRDS(keratinocytes, file = rds_kc)
saveRDS(krt35_85_kc,   file = rds_krt35_85)
saveRDS(krt35_85_DEGs, file = rds_krt35_DEG)
saveRDS(krt31_clusters, file = file.path(kc_dir, "krt31_clusters.rds"))

cat(sprintf("  Saved: %s\n", rds_kc))
cat(sprintf("  Saved: %s\n", rds_krt35_85))
cat(sprintf("  Saved: %s\n", rds_krt35_DEG))

cat("\n=== Keratinocyte analysis complete ===\n")
