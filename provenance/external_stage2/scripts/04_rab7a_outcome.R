# Stage 2 (2026-10-01): RAB7A outcome in the frozen population (plan section 6).
# GATED: run only after checkpoint C (population freeze) has been independently
# reviewed. Pass --reviewed to confirm. Run from the project root.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_common.R"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(dplyr) })
if (!"--reviewed" %in% commandArgs(trailingOnly = TRUE))
  stop("Checkpoint C review not confirmed; rerun with --reviewed after independent review")
plan_sha <- check_plan_frozen()
fz <- fromJSON(file.path(out_dir, "C9_population_freeze.json"))
stopifnot(identical(fz$plan_sha256, plan_sha))
if (fz$status == "no_cluster_qualifies_inconclusive") stop("No population (plan 5.2): inconclusive, not run")
kc_rds <- file.path(derived, "stage2_keratinocytes_subclustered.rds")
stopifnot(identical(digest(kc_rds, algo = "sha256", file = TRUE), fz$kc_object_sha256))
cells <- readLines(file.path(derived, "stage2_population_cells.txt"))
stopifnot(identical(digest(paste(cells, collapse = "\n"), algo = "sha256", serialize = FALSE),
                    fz$population_cells_sha256))
out <- file.path(out_dir, "D_rab7a")
if (dir.exists(out)) stop("Outcome already produced; refusing to overwrite ", out)
dir.create(out)

kc <- readRDS(kc_rds)
pop <- subset(kc, cells = cells); rm(kc); invisible(gc())
cnt <- LayerData(pop, assay = "RNA", layer = "counts")
dat <- LayerData(pop, assay = "RNA", layer = "data")
md <- pop@meta.data
md$RAB7A_counts <- cnt["RAB7A", ]; md$total_counts <- colSums(cnt)
md$RAB7A_lognorm <- dat["RAB7A", ]
rm(cnt, dat)

# Primary + co-primary, per sample
ps <- md |> group_by(sample, phenotype, donor) |>
  summarise(n_cells = n(), rab7a_umi = sum(RAB7A_counts), total_umi = sum(total_counts),
            pseudobulk_log2cpm = log2(sum(RAB7A_counts) / sum(total_counts) * 1e6 + 1),
            mean_lognorm = mean(RAB7A_lognorm), detection = mean(RAB7A_counts > 0),
            .groups = "drop")
write.table(ps, file.path(out, "D1_rab7a_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Primary: per paired donor, White - Black (no P value; plan 6)
pair <- do.call(rbind, lapply(c("F31", "F62"), function(d) {
  b <- ps[ps$donor == d & ps$phenotype == "Black", ]; w <- ps[ps$donor == d & ps$phenotype == "White", ]
  data.frame(donor = d,
             diff_pseudobulk_log2cpm_white_minus_black = w$pseudobulk_log2cpm - b$pseudobulk_log2cpm,
             diff_mean_lognorm_white_minus_black = w$mean_lognorm - b$mean_lognorm,
             n_black = b$n_cells, n_white = w$n_cells)
}))
dirs <- sign(pair$diff_pseudobulk_log2cpm_white_minus_black)
verdict <- if (all(dirs < 0)) "Direction consistent in both paired donors (white < black)" else
  if (all(dirs > 0)) "Opposite in both paired donors (white > black)" else "Inconsistent across donors"
write.table(pair, file.path(out, "D2_rab7a_paired_donors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Secondary: sample-level exact Wilcoxon, 4 Black vs 2 White (min two-sided p = 2/15)
wt <- wilcox.test(pseudobulk_log2cpm ~ phenotype, data = ps, exact = TRUE)
# Descriptive only: cell-level Wilcoxon (pseudoreplicated)
wc <- wilcox.test(RAB7A_lognorm ~ phenotype, data = md, exact = FALSE)

pdf(file.path(out, "D3_rab7a_violin_by_sample_DESCRIPTIVE.pdf"), width = 9, height = 5)
print(VlnPlot(pop, "RAB7A", group.by = "sample", pt.size = 0) +
        ggplot2::labs(subtitle = "Descriptive only; frozen population; LogNormalize of genuine counts"))
dev.off()

write_json(list(script = "04_rab7a_outcome.R", plan_sha256 = plan_sha, population_status = fz$status,
                primary_verdict = verdict, paired = pair,
                secondary_sample_level_wilcoxon = list(W = unname(wt$statistic), p = wt$p.value,
                  note = "exact, 4 Black vs 2 White samples; minimum attainable two-sided p = 2/15"),
                descriptive_cell_level_wilcoxon = list(W = unname(wc$statistic), p = wc$p.value,
                  n_cells = nrow(md), note = "pseudoreplicated; not used for inference"),
                session = capture.output(sessionInfo())),
           file.path(out, "D0_rab7a_outcome.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(ps); print(pair); cat(verdict, "\n")
