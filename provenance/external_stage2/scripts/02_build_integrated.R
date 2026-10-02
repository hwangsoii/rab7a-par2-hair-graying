# Stage 2 (2026-10-01): QC, CCA integration, clustering and the deterministic
# marker-based annotation (plan sections 3-4). Input: Zenodo 15103193
# Cell Ranger-filtered matrices only. Does not read or report RAB7A/F2RL1.
# Run from the project root after checkpoint A passed and the plan was frozen.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
plan_sha <- check_plan_frozen()
set.seed(1337)
# Technical only (not an analysis parameter): allow FindIntegrationAnchors to
# pass the ~0.9 GB object list to future's (sequential) workers.
options(future.globals.maxSize = 8 * 1024^3)

# ---- 3.1-3.3 per-sample load, QC, normalisation (sequential, sparse) --------
qc <- list(); obj_list <- list()
for (i in seq_len(nrow(samples))) {
  s <- samples$sample[i]
  p <- function(x) file.path(zen_dir, sprintf("%s_filtered_%s", s, x))
  m <- as(readMM(gzfile(p("matrix.mtx.gz"))), "CsparseMatrix")
  feats <- read.delim(gzfile(p("features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  rownames(m) <- make.unique(feats$V2)
  colnames(m) <- paste0(s, "_", readLines(gzfile(p("barcodes.tsv.gz"))))
  n_cr <- ncol(m)
  x <- CreateSeuratObject(counts = m, project = s, min.cells = 3, min.features = 200)
  rm(m); invisible(gc())
  x$sample <- s; x$phenotype <- samples$phenotype[i]; x$donor <- samples$donor[i]
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  n0 <- ncol(x)
  md <- x@meta.data
  fail <- data.frame(nFeature_le_500 = sum(md$nFeature_RNA <= 500),
                     nFeature_ge_5000 = sum(md$nFeature_RNA >= 5000),
                     percent_mt_ge_10 = sum(md$percent.mt >= 10),
                     nCount_ge_40000 = sum(md$nCount_RNA >= 40000))
  x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 &
                percent.mt < 10 & nCount_RNA < 40000)
  qc[[s]] <- data.frame(sample = s, cells_cellranger_filtered = n_cr,
                        cells_min_features_200 = n0, fail,
                        cells_after_qc = ncol(x))
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj_list[[s]] <- x
  cat(sprintf("%s: %d -> %d cells\n", s, n0, ncol(x)))
}
qc <- do.call(rbind, qc)
write.table(qc, file.path(out_dir, "C1_qc_cells_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 3.4-3.5 CCA integration and clustering -------------------------------
anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:20, verbose = FALSE)
obj <- IntegrateData(anchorset = anchors, dims = 1:20, verbose = FALSE)
rm(anchors, obj_list); invisible(gc())
DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:20, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)

# ---- 4. deterministic marker-based annotation -------------------------------
lineages <- list(
  KC_basal_ORS = c("KRT5", "KRT14", "KRT15", "KRT17", "KRT6A", "KRT16", "COL17A1"),
  KC_IFE_suprabasal = c("KRT1", "KRT10", "KRTDAP", "SBSN", "DMKN"),
  KC_IRS = c("KRT25", "KRT27", "KRT71", "KRT73", "TCHH"),
  KC_cortex_cuticle = c("KRT31", "KRT35", "KRT85", "KRT32", "KRT33A"),
  KC_matrix = c("MSX2", "LEF1", "HOXC13"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  Lymphatic = c("LYVE1", "PROX1", "CCL21"),
  Mural = c("RGS5", "ACTA2", "TAGLN", "MYH11", "PDGFRB"),
  Melanocyte = c("PMEL", "MLANA", "DCT", "TYRP1", "TYR"),
  T_NK = c("PTPRC", "CD3D", "CD3E", "CD2", "NKG7"),
  Myeloid = c("CD68", "CD14", "LYZ", "AIF1", "CD207"),
  Mast = c("TPSAB1", "CPA3", "MS4A2"),
  Gland = c("DCD", "SCGB2A2", "PIP", "MUCL1"),
  Schwann = c("MPZ", "PLP1", "PRX"))
stopifnot(!any(outcome_genes %in% unlist(lineages)))
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
missing <- lapply(lineages, function(g) setdiff(g, rownames(obj)))
lineages <- lapply(lineages, function(g) intersect(g, rownames(obj)))
obj <- AddModuleScore(obj, features = lineages, name = "lin_", assay = "RNA", seed = 1)
sc_cols <- paste0("lin_", seq_along(lineages))
cl_scores <- obj@meta.data |> group_by(seurat_clusters) |>
  summarise(across(all_of(sc_cols), mean), n_cells = n(), .groups = "drop")
mat <- as.matrix(cl_scores[, sc_cols]); colnames(mat) <- names(lineages)
ord <- t(apply(mat, 1, order, decreasing = TRUE))
annot <- data.frame(cluster = cl_scores$seurat_clusters, n_cells = cl_scores$n_cells,
                    top_lineage = names(lineages)[ord[, 1]],
                    second_lineage = names(lineages)[ord[, 2]],
                    margin = round(mat[cbind(seq_len(nrow(mat)), ord[, 1])] -
                                   mat[cbind(seq_len(nrow(mat)), ord[, 2])], 4),
                    round(mat, 4), check.names = FALSE)
annot$keratinocyte <- grepl("^KC_", annot$top_lineage)
obj$lineage <- annot$top_lineage[match(obj$seurat_clusters, annot$cluster)]
obj$keratinocyte <- grepl("^KC_", obj$lineage)
write.table(annot, file.path(out_dir, "C2_cluster_annotation_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
comp <- as.data.frame.matrix(table(obj$seurat_clusters, obj$sample))
write.table(cbind(cluster = rownames(comp), comp), file.path(out_dir, "C3_cluster_by_sample_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(obj, file.path(derived, "stage2_integrated_annotated.rds"))
pdf(file.path(out_dir, "C4_umap_clusters_lineage_sample.pdf"), width = 21, height = 6.5)
print(DimPlot(obj, group.by = "seurat_clusters", label = TRUE) + NoLegend() |
      DimPlot(obj, group.by = "lineage", label = TRUE, repel = TRUE) |
      DimPlot(obj, group.by = "sample"))
dev.off()

write_json(list(script = "02_build_integrated.R", plan_sha256 = plan_sha,
                n_cells = ncol(obj), n_clusters = nlevels(obj$seurat_clusters),
                keratinocyte_cells = sum(obj$keratinocyte),
                keratinocyte_clusters = as.character(annot$cluster[annot$keratinocyte]),
                missing_marker_genes = missing,
                session = capture.output(sessionInfo())),
           file.path(out_dir, "C0_build_summary.json"), auto_unbox = TRUE, pretty = TRUE)
print(annot[, c("cluster", "n_cells", "top_lineage", "second_lineage", "margin")])
