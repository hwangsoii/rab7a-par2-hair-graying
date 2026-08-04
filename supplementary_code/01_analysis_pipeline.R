# =============================================================================
# 01_analysis_pipeline.R — QC, Integration, Clustering, Annotation
#
# Corrected version of Script_Totalcell_hair_ver1.1_swh.Rmd
# Fixes:
#   1. QC filtering applied to ALL 4 samples (bug: black2/grey2 were unfiltered)
#   2. Correct remove() calls (was using wrong variable names)
#   3. Uses config.R for paths (no more hard-coded absolute paths)
#
# Input:  10x Cell Ranger filtered_feature_bc_matrix (4 samples)
# Output: data.hair.std.r0.5.rds, data.hair.rename.rds
# =============================================================================

source("config.R")

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# =============================================================================
# STEP 1: Load raw data
# =============================================================================
cat(">>> Loading 10x data ...\n")

require_sample_dirs(c("Black1", "Black2", "Grey1", "Grey2"))

black1.data <- Read10X(data.dir = file.path(raw_dir, "Black1/filtered_feature_bc_matrix"))
black2.data <- Read10X(data.dir = file.path(raw_dir, "Black2/filtered_feature_bc_matrix"))
grey1.data  <- Read10X(data.dir = file.path(raw_dir, "Grey1/filtered_feature_bc_matrix"))
grey2.data  <- Read10X(data.dir = file.path(raw_dir, "Grey2/filtered_feature_bc_matrix"))

black1 <- CreateSeuratObject(counts = black1.data, project = "black",  min.cells = 3, min.features = 200)
black2 <- CreateSeuratObject(counts = black2.data, project = "black2", min.cells = 3, min.features = 200)
grey1  <- CreateSeuratObject(counts = grey1.data,  project = "grey",   min.cells = 3, min.features = 200)
grey2  <- CreateSeuratObject(counts = grey2.data,  project = "grey2",  min.cells = 3, min.features = 200)

rm(black1.data, black2.data, grey1.data, grey2.data)

# =============================================================================
# STEP 2: QC — percent.mt and filtering (ALL 4 samples)
# =============================================================================
cat(">>> QC filtering ...\n")

samples <- list(black1 = black1, black2 = black2, grey1 = grey1, grey2 = grey2)

for (nm in names(samples)) {
  samples[[nm]][["percent.mt"]] <- PercentageFeatureSet(samples[[nm]], pattern = "^MT-")
  n_before <- ncol(samples[[nm]])
  samples[[nm]] <- subset(samples[[nm]],
                          subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
                                   percent.mt < 10 & nCount_RNA < 40000)
  n_after <- ncol(samples[[nm]])
  cat(sprintf("  %s: %d -> %d cells (removed %d)\n", nm, n_before, n_after, n_before - n_after))
}

black1 <- samples$black1
black2 <- samples$black2
grey1  <- samples$grey1
grey2  <- samples$grey2
rm(samples)

# =============================================================================
# STEP 3: Protocol labels, normalization, feature selection
# =============================================================================
cat(">>> Normalizing and finding variable features ...\n")

black1$protocol <- "Black"
black2$protocol <- "Black"
grey1$protocol  <- "Grey"
grey2$protocol  <- "Grey"

obj_list <- list(black1, grey1, black2, grey2)
obj_list <- lapply(obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x
})

# =============================================================================
# STEP 4: CCA Integration
# =============================================================================
cat(">>> Finding integration anchors (dims 1:20) ...\n")
anchors  <- FindIntegrationAnchors(object.list = obj_list, dims = 1:20)
hair.std <- IntegrateData(anchorset = anchors, dims = 1:20)
rm(obj_list, anchors, black1, black2, grey1, grey2); gc()

# =============================================================================
# STEP 5: Scaling, PCA, UMAP, Clustering
# =============================================================================
cat(">>> Scaling, PCA, UMAP, clustering ...\n")

DefaultAssay(hair.std) <- "integrated"
hair.std <- ScaleData(hair.std)
hair.std <- RunPCA(hair.std, npcs = 30)
hair.std <- RunUMAP(hair.std, reduction = "pca", dims = 1:20)
hair.std <- RunTSNE(hair.std, reduction = "pca", dims = 1:20)
hair.std <- FindNeighbors(hair.std, reduction = "pca", dims = 1:20)

hair.std.r0.5 <- FindClusters(hair.std, resolution = 0.5)
rm(hair.std)

cat(">>> Saving clustered object ...\n")
saveRDS(hair.std.r0.5, file = rds_std)

# =============================================================================
# STEP 6: Cell type annotation
# =============================================================================
cat(">>> Annotating cell types ...\n")

rename <- hair.std.r0.5
rename <- RenameIdents(rename,
                       `0`  = "Fibroblast",
                       `1`  = "Lymphocyte",
                       `2`  = "Keratinocyte",
                       `3`  = "Keratinocyte",
                       `4`  = "Vascular",
                       `5`  = "Keratinocyte",
                       `6`  = "Vascular",
                       `7`  = "Myeloid",
                       `8`  = "Keratinocyte",
                       `9`  = "Pericyte",
                       `10` = "Keratinocyte",
                       `11` = "Keratinocyte",
                       `12` = "Keratinocyte",
                       `13` = "Melanocyte",
                       `14` = "Pericyte",
                       `15` = "Myeloid",
                       `16` = "Fibroblast",
                       `17` = "Mast cells/Basophils",
                       `18` = "Gland",
                       `19` = "Lymphocyte")

saveRDS(rename, file = rds_rename)

cat(">>> Cell type counts:\n")
print(table(Idents(rename)))

cat("\n=== Pipeline complete ===\n")
cat(sprintf("  Saved: %s\n", rds_std))
cat(sprintf("  Saved: %s\n", rds_rename))
