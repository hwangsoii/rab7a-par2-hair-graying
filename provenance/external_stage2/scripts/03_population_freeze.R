# Stage 2 (2026-10-01): keratinocyte subclustering and the single frozen
# population rule (plan section 5). Never reads RAB7A or F2RL1.
# Output: checkpoint C freeze files. Run from the project root.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
plan_sha <- check_plan_frozen()
freeze_json <- file.path(out_dir, "C9_population_freeze.json")
if (file.exists(freeze_json)) stop("Population already frozen; refusing to overwrite ", freeze_json)

obj <- readRDS(file.path(derived, "stage2_integrated_annotated.rds"))
kc <- subset(obj, subset = keratinocyte)
rm(obj); invisible(gc())

# 5.1 mirrors 02_keratinocyte_analysis.R STEP 1
DefaultAssay(kc) <- "RNA"
kc[["integrated"]] <- NULL
kc <- NormalizeData(kc, verbose = FALSE)
kc <- FindVariableFeatures(kc, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
kc <- ScaleData(kc, verbose = FALSE)
kc <- RunPCA(kc, npcs = 30, verbose = FALSE)
kc <- FindNeighbors(kc, reduction = "pca", dims = 1:20, verbose = FALSE)
kc <- FindClusters(kc, resolution = 0.6, verbose = FALSE)
kc <- RunUMAP(kc, reduction = "pca", dims = 1:20, n.neighbors = 30, min.dist = 0.1,
              spread = 1.2, verbose = FALSE)

# 5.2 the single population rule: cluster mean of (KRT35 + KRT85)/2 > 0.5
dat <- LayerData(kc, assay = "RNA", layer = "data")
kc$KRT35_KRT85_score <- (dat["KRT35", ] + dat["KRT85", ]) / 2
rm(dat)
cl <- kc@meta.data |> group_by(kc_cluster = seurat_clusters) |>
  summarise(mean_krt35_85 = mean(KRT35_KRT85_score), n_cells = n(), .groups = "drop") |>
  arrange(desc(mean_krt35_85))
sel <- as.character(cl$kc_cluster[cl$mean_krt35_85 > 0.5])
status <- if (length(sel) == 0) "no_cluster_qualifies_inconclusive" else
  if (length(sel) == nrow(cl)) "rule_non_discriminating_all_keratinocytes" else "selected_subset"
kc$population <- as.character(kc$seurat_clusters) %in% sel
cl$selected <- as.character(cl$kc_cluster) %in% sel
by_sample <- as.data.frame.matrix(table(kc$seurat_clusters, kc$sample))
write.table(cbind(cl, by_sample[as.character(cl$kc_cluster), samples$sample]),
            file.path(out_dir, "C5_kc_subcluster_krt35_85_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
lin <- as.data.frame.matrix(table(kc$seurat_clusters, kc$lineage))
write.table(cbind(kc_cluster = rownames(lin), lin), file.path(out_dir, "C6_kc_subcluster_by_lineage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
pop <- kc@meta.data |> group_by(sample, phenotype, donor) |>
  summarise(keratinocytes = n(), population_cells = sum(population), .groups = "drop")
write.table(pop, file.path(out_dir, "C7_population_cells_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(out_dir, "C8_kc_subclusters_population.pdf"), width = 21, height = 6.5)
print(DimPlot(kc, group.by = "seurat_clusters", label = TRUE) + NoLegend() |
      FeaturePlot(kc, "KRT35_KRT85_score") |
      DimPlot(kc, group.by = "population") | DimPlot(kc, group.by = "sample"))
dev.off()

cells <- sort(colnames(kc)[kc$population])
writeLines(cells, file.path(derived, "stage2_population_cells.txt"))
saveRDS(kc, file.path(derived, "stage2_keratinocytes_subclustered.rds"))
write_json(list(script = "03_population_freeze.R", frozen_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                plan_sha256 = plan_sha, rule = "KC subcluster mean (KRT35+KRT85)/2 > 0.5 (RNA data layer)",
                status = status, n_kc = ncol(kc), n_kc_clusters = nrow(cl),
                selected_clusters = sel, n_population_cells = length(cells),
                population_cells_sha256 = digest(paste(cells, collapse = "\n"), algo = "sha256", serialize = FALSE),
                outcome_script_sha256 = digest(file.path(stage2_dir, "scripts", "04_rab7a_outcome.R"), algo = "sha256", file = TRUE),
                kc_object_sha256 =digest(file.path(derived, "stage2_keratinocytes_subclustered.rds"), algo = "sha256", file = TRUE),
                session = capture.output(sessionInfo())),
           freeze_json, auto_unbox = TRUE, pretty = TRUE)
print(as.data.frame(cl)); print(as.data.frame(pop)); cat("status:", status, "\n")
