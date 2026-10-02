# Amendment 01 (2026-10-01): marker-only annotation-consistency review of the
# v2 candidate (existing lineage == "KC_cortex_cuticle" in the integrated
# object; amendment section 3). Extracts only the predefined review markers;
# never reads RAB7A or F2RL1. Writes aggregate tables/plots (outputs/V2_M*)
# and no per-cell data or object. Run from the project root.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_amend_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2); library(patchwork) })
plan_sha <- check_plan_frozen()
amend_sha <- check_amendment_frozen()
a4_sha <- check_a4_pass()
if (length(list.files(out_dir, pattern = "^V2_M"))) stop("V2 marker outputs already exist; refusing to overwrite")
log <- character()
note <- function(...) { m <- paste0(...); log <<- c(log, m); cat(m, "\n") }

obj_sha <- sha_file(integrated_rds)
note("integrated object sha256 before load: ", obj_sha)
obj <- readRDS(integrated_rds)
md <- obj@meta.data
stopifnot(all(c("lineage", "seurat_clusters", "sample", "keratinocyte") %in% colnames(md)))
md$cluster <- as.character(md$seurat_clusters)

# ---- S1 membership: derived from the existing label --------------------------
c2 <- read.delim(file.path(out_dir, "C2_cluster_annotation_scores.tsv"), stringsAsFactors = FALSE)
lab_by_cl <- tapply(md$lineage, md$cluster, function(x) paste(unique(x), collapse = "|"))
recon <- data.frame(cluster = as.character(c2$cluster), c2_top_lineage = c2$top_lineage,
                    object_lineage = unname(lab_by_cl[as.character(c2$cluster)]),
                    c2_n = c2$n_cells, object_n = as.integer(table(md$cluster)[as.character(c2$cluster)]))
recon$agree <- recon$c2_top_lineage == recon$object_lineage & recon$c2_n == recon$object_n
note("C2 reconciliation: ", sum(recon$agree), "/", nrow(recon), " clusters agree")
stopifnot(all(recon$agree))

md$candidate <- md$lineage == v2_label
cand_clusters <- sort(as.integer(unique(md$cluster[md$candidate])))
cand_by_sample <- table(factor(md$sample[md$candidate], levels = samples$sample))
audit <- c(F18 = 179, F59 = 820, F31B = 547, F31W = 247, F62B = 1138, F62W = 384)
s1 <- data.frame(sample = samples$sample, derived = as.integer(cand_by_sample[samples$sample]),
                 audit_value = unname(audit[samples$sample]))
s1$matches_audit <- s1$derived == s1$audit_value
note("candidate clusters (integrated res 0.5): ", paste(cand_clusters, collapse = ","),
     "; n = ", sum(md$candidate), "; audit clusters 3,6,9,14 / n 3315 match: ",
     identical(cand_clusters, c(3L, 6L, 9L, 14L)) && sum(md$candidate) == 3315)
v1_cells <- readLines(v1_cells_txt)
c9 <- fromJSON(c9_json)
stopifnot(identical(digest(paste(v1_cells, collapse = "\n"), algo = "sha256", serialize = FALSE),
                    c9$population_cells_sha256))
cand_cells <- rownames(md)[md$candidate]
note("candidate subset of v1 (C9 list): ", all(cand_cells %in% v1_cells))
stopifnot(all(cand_cells %in% v1_cells), all(md$keratinocyte[md$candidate]))

# ---- marker extraction (predefined markers only) ----------------------------
genes_req <- unlist(review_markers, use.names = FALSE)
set_of <- setNames(rep(names(review_markers), lengths(review_markers)), genes_req)
genes <- intersect(genes_req, rownames(obj[["RNA"]]))
note("markers absent from RNA assay: ", length(setdiff(genes_req, genes)), " (", paste(setdiff(genes_req, genes), collapse = ","), ")")
stopifnot(!any(outcome_genes %in% genes))
kc_cells <- rownames(md)[md$keratinocyte]
dat <- LayerData(obj, assay = "RNA", layer = "data")[genes, kc_cells, drop = FALSE]
cnt <- LayerData(obj, assay = "RNA", layer = "counts")[genes, kc_cells, drop = FALSE]
mk <- md[kc_cells, ]
det <- cnt > 0
rm(cnt); invisible(gc())

summarise_groups <- function(groups, keys) {
  idx <- split(seq_along(groups), groups, drop = TRUE)
  do.call(rbind, lapply(names(idx), function(g) {
    i <- idx[[g]]
    data.frame(keys[i[1], , drop = FALSE], n_cells = length(i), gene = genes, marker_set = set_of[genes],
               mean_lognorm = round(Matrix::rowMeans(dat[, i, drop = FALSE]), 4),
               detection = round(Matrix::rowMeans(det[, i, drop = FALSE]), 4), row.names = NULL)
  }))
}
w <- function(x, f) write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)

m1 <- summarise_groups(mk$cluster, data.frame(cluster = mk$cluster, lineage = mk$lineage, candidate = mk$candidate))
m1 <- m1[order(!m1$candidate, as.integer(m1$cluster)), ]
w(m1, "V2_M1_markers_by_integrated_cluster.tsv")
m2 <- summarise_groups(mk$lineage, data.frame(lineage = mk$lineage))
w(m2, "V2_M2_markers_by_lineage_group.tsv")
grp3 <- paste(mk$sample, ifelse(mk$candidate, "candidate", "other_keratinocytes"))
m3 <- summarise_groups(grp3, data.frame(sample = mk$sample, group = ifelse(mk$candidate, "candidate", "other_keratinocytes")))
w(m3, "V2_M3_markers_candidate_vs_other_by_sample.tsv")
cand_grp <- function(x) ifelse(mk$candidate, x, NA)  # full-length groups; NA = not candidate
m4 <- summarise_groups(cand_grp(paste(mk$cluster, mk$sample)),
                       data.frame(cluster = mk$cluster, sample = mk$sample))
w(m4, "V2_M4_markers_candidate_cluster_by_sample.tsv")

# S3 co-detection profiles (fractions of cells), by keratinocyte cluster and,
# for the candidate, by sample
g <- function(x) as.logical(det[x, ])
hk <- g("KRT35") | g("KRT85")
prof <- data.frame(KRT31_pos = g("KRT31"), KRT35_pos_KRT31_neg = g("KRT35") & !g("KRT31"),
                   KRT35or85_pos = hk, KRT35or85_and_KRT14 = hk & g("KRT14"),
                   KRT35or85_and_GATA3_TCHH = hk & g("GATA3") & g("TCHH"),
                   KRT35or85_and_MSX2orLEF1 = hk & (g("MSX2") | g("LEF1")),
                   GATA3_and_TCHH = g("GATA3") & g("TCHH"), MSX2_or_LEF1 = g("MSX2") | g("LEF1"),
                   MKI67_pos = g("MKI67"), KRT14_pos = g("KRT14"))
codet <- function(groups, keys) do.call(rbind, lapply(split(seq_along(groups), groups, drop = TRUE), function(i)
  data.frame(keys[i[1], , drop = FALSE], n_cells = length(i), t(round(colMeans(prof[i, ]), 4)), row.names = NULL)))
m5a <- codet(mk$cluster, data.frame(level = "cluster", cluster = mk$cluster, sample = "all",
                                     lineage = mk$lineage, candidate = mk$candidate))
m5b <- codet(cand_grp(paste(mk$cluster, mk$sample)),
             data.frame(level = "candidate_cluster_x_sample", cluster = mk$cluster, sample = mk$sample,
                        lineage = mk$lineage, candidate = mk$candidate))
m5c <- codet(cand_grp(mk$sample), data.frame(level = "candidate_pooled_by_sample", cluster = "candidate",
                                             sample = mk$sample, lineage = mk$lineage, candidate = mk$candidate))
m5 <- rbind(m5a[order(!m5a$candidate, as.integer(m5a$cluster)), ], m5c, m5b)
w(m5, "V2_M5_codetection_profiles.tsv")
s1_out <- list(derived_clusters = cand_clusters, n = sum(md$candidate), per_sample = s1,
               all_six_samples_present = all(s1$derived > 0), c2_reconciliation = recon)
w(s1, "V2_M7_S1_membership_vs_audit.tsv")

# ---- dot plots (RNA data layer; display only, not evidence of identity) -------
kc <- subset(obj, cells = kc_cells); rm(obj); invisible(gc())
DefaultAssay(kc) <- "RNA"
kc$cluster_lab <- factor(paste0(kc$seurat_clusters, ifelse(kc$lineage == v2_label, "* ", " "), kc$lineage),
                         levels = unique(paste0(m1$cluster, ifelse(m1$candidate, "* ", " "), m1$lineage)))
feat <- lapply(review_markers, intersect, genes)
p1 <- DotPlot(kc, features = feat, group.by = "cluster_lab", assay = "RNA") + RotatedAxis() +
  labs(title = "Review markers by original integrated cluster (keratinocytes; * = KC_cortex_cuticle candidate)",
       subtitle = "Descriptive only. Colour = scaled mean LogNormalize; size = detection fraction") +
  theme(axis.text.x = element_text(size = 7))
cand <- subset(kc, cells = cand_cells)
cand$cl_sample <- factor(paste0("cl", cand$seurat_clusters, " ", cand$sample),
                         levels = as.vector(outer(paste0("cl", cand_clusters), samples$sample, paste)))
cand$sample_f <- factor(cand$sample, levels = samples$sample)
p2 <- DotPlot(cand, features = feat, group.by = "sample_f", assay = "RNA") + RotatedAxis() +
  labs(title = "Candidate, pooled, by sample") + theme(axis.text.x = element_text(size = 7))
p3 <- DotPlot(cand, features = feat, group.by = "cl_sample", assay = "RNA") + RotatedAxis() +
  labs(title = "Candidate cluster x sample") + theme(axis.text.x = element_text(size = 7))
pdf(file.path(out_dir, "V2_M6_marker_dotplots.pdf"), width = 16, height = 8)
print(p1); print(p2); print(p3)
dev.off()

outs <- list.files(out_dir, pattern = "^V2_M[1-7]_", full.names = TRUE)
write_json(list(script = "06_marker_review_v2.R", run_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                plan_sha256 = plan_sha, amendment_sha256 = amend_sha, a4_sha256 = a4_sha,
                integrated_object_sha256 = obj_sha,
                integrated_object_sha256_after = sha_file(integrated_rds),
                script06_sha256 = sha_file(file.path(stage2_dir, "scripts", "06_marker_review_v2.R")),
                markers_requested = review_markers, markers_absent = setdiff(genes_req, genes),
                outcome_genes_extracted = FALSE, s1 = s1_out, log = log,
                output_sha256 = as.list(setNames(vapply(outs, sha_file, ""), basename(outs))),
                session = capture.output(sessionInfo())),
           file.path(out_dir, "V2_M0_marker_review.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
