# Amendment 01 (2026-10-01): descriptive RAB7A outcome in the frozen v1
# (contextual, all keratinocytes) and v2 (targeted, KC_cortex_cuticle)
# populations (amendment section 5). GATED (amendment section 6):
#   1. --review-record=<file> must exist, approved == true, >= 2 reviewers, ISO date
#   2. current C10 SHA-256 == record$c10_sha256 (C10 holds no self-checksum)
#   3. every C10 dependency matches the current file AND record$approved_versions
#   4. A4 status PASS
# --gate-only stops after the gate, before readRDS, gene extraction or any
# output-directory creation. Run from the project root.
# Usage: Rscript 08_outcome_amendment01.R --review-record=<file> [--gate-only]
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_amend_common.R"))
args <- commandArgs(trailingOnly = TRUE)
gate_only <- "--gate-only" %in% args
stage <- function(s) cat("GATE_STAGE:", s, "\n")
gate_fail <- function(reason, detail = "") {
  cat("GATE_FAIL:", reason, if (nzchar(detail)) paste0("(", detail, ")") else "", "\n")
  quit(save = "no", status = 2)
}
required_deps <- c("plan", "B0", "C9", "amendment", "B1", "v1_cells", "v2_cells", "kc_object",
                   "integrated_object", "A4", "V2_M0", "V2_D0", "C2_report", "script05", "script06",
                   "script07", "script08", "stage2_common", "stage2_validate", "stage2_amend_common")
out <- file.path(out_dir, "E_outcome_amendment01")

# ---- gate ------------------------------------------------------------------
stage("approval_record")
rr_path <- sub("^--review-record=", "", grep("^--review-record=", args, value = TRUE))
if (length(rr_path) != 1 || !file.exists(rr_path)) gate_fail("review_record_missing")
rr <- tryCatch(fromJSON(rr_path, simplifyVector = FALSE), error = function(e) NULL)
if (is.null(rr)) gate_fail("review_record_unreadable")
if (!isTRUE(rr$approved) || length(rr$reviewers) < 2 || !is.character(rr$review_date) ||
    !grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}", rr$review_date)) gate_fail("review_not_approved")
stage("c10")
if (!file.exists(c10_json)) gate_fail("c10_missing")
if (!identical(sha_file(c10_json), rr$c10_sha256)) gate_fail("c10_mismatch")
c10 <- fromJSON(c10_json, simplifyVector = FALSE)
deps <- c10$dependencies
dep_names <- vapply(deps, function(d) as.character(d$name), "")
if (length(setdiff(required_deps, dep_names))) gate_fail("c10_missing_dependency", paste(setdiff(required_deps, dep_names), collapse = ","))
stage("dependencies")
v <- validate_dependencies(deps, base_dir = ".")
if (!v$ok) gate_fail("dependency_mismatch", paste(v$failures$name, v$failures$reason, collapse = ";"))
for (d in deps) if (!identical(rr$approved_versions[[d$name]], d$sha256))
  gate_fail(paste0("dependency_not_approved:", d$name))
if (!identical(fromJSON(a4_json)$status, "PASS")) gate_fail("a4_not_pass")
if (dir.exists(out)) gate_fail("outcome_folder_exists")
stage("gate_ok"); cat("GATE_OK\n")
if (gate_only) quit(save = "no", status = 0)

# ---- outcome (reached only with a genuine approval record) ------------------
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })
category <- function(x) if (!is.finite(x)) stop("non-finite difference") else
  if (x < 0) "white lower" else if (x > 0) "white higher" else "zero difference"

summarise_population <- function(label, rds, cells) {
  obj <- readRDS(rds)
  if (!all(cells %in% colnames(obj))) stop(label, ": cells missing from source object")
  md <- obj@meta.data[cells, ]
  cnt <- LayerData(obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
  if (!all(cnt@x >= 0 & cnt@x == round(cnt@x))) stop(label, ": counts not nonnegative integers")
  lib <- Matrix::colSums(cnt)
  if (!isTRUE(all.equal(unname(lib), md$nCount_RNA, tolerance = 0))) stop(label, ": colSums(counts) != nCount_RNA")
  feature_universe <- nrow(cnt)
  rab_cnt <- cnt["RAB7A", ]; rm(cnt)
  rab_dat <- LayerData(obj, assay = "RNA", layer = "data")["RAB7A", cells]
  rm(obj); invisible(gc())
  got <- sort(unique(md$sample))
  if (!identical(got, sort(samples$sample))) stop(label, ": sample set ", paste(got, collapse = ","), " != intended six")
  ps <- do.call(rbind, lapply(samples$sample, function(s) {
    i <- md$sample == s
    den <- sum(lib[i]); num <- sum(rab_cnt[i])
    if (sum(i) == 0 || !is.finite(den) || den <= 0) stop(label, ": sample ", s, " has 0 cells or invalid denominator")
    data.frame(population = label, sample = s, phenotype = samples$phenotype[samples$sample == s],
               donor = samples$donor[samples$sample == s], n_cells = sum(i), rab7a_umi = num,
               denominator_umi = den, pseudobulk_log2cpm1 = log2(num / den * 1e6 + 1),
               mean_lognorm = mean(rab_dat[i]), detection = mean(rab_cnt[i] > 0))
  }))
  pair <- do.call(rbind, lapply(c("F31", "F62"), function(d) {
    b <- ps[ps$donor == d & ps$phenotype == "Black", ]; w <- ps[ps$donor == d & ps$phenotype == "White", ]
    if (nrow(b) != 1 || nrow(w) != 1) stop(label, ": donor ", d, " lacks exactly one Black and one White sample")
    dp <- w$pseudobulk_log2cpm1 - b$pseudobulk_log2cpm1; dm <- w$mean_lognorm - b$mean_lognorm
    data.frame(population = label, donor = d, n_black = b$n_cells, n_white = w$n_cells,
               diff_pseudobulk_log2cpm1 = dp, category_pseudobulk = category(dp),
               diff_mean_lognorm = dm, category_mean_lognorm = category(dm),
               metric_disagreement = sign(dp) != sign(dm))
  }))
  cells_df <- data.frame(population = label, sample = md$sample, rab7a_lognorm = rab_dat,
                         detected = rab_cnt > 0)
  list(ps = ps, pair = pair, cells = cells_df, feature_universe = feature_universe)
}

v1 <- summarise_population("v1_all_keratinocytes", kc_rds, readLines(v1_cells_txt))
v2 <- summarise_population("v2_KC_cortex_cuticle", integrated_rds, readLines(v2_cells_txt))
if (!dir.create(out, showWarnings = FALSE)) stop("Refusing to overwrite ", out)
ps <- rbind(v1$ps, v2$ps); pair <- rbind(v1$pair, v2$pair)
write.table(ps, file.path(out, "E1_per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pair, file.path(out, "E2_paired_donors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
all_lower <- tapply(pair$category_pseudobulk == "white lower" & pair$category_mean_lognorm == "white lower",
                    pair$population, all)
cells <- rbind(v1$cells, v2$cells); cells$sample <- factor(cells$sample, levels = samples$sample)
pdf(file.path(out, "E3_distribution_detection_DESCRIPTIVE.pdf"), width = 11, height = 7.5)
print(ggplot(cells, aes(sample, rab7a_lognorm)) + geom_violin(scale = "width") + facet_wrap(~population, ncol = 1) +
        theme_classic() + labs(y = "RAB7A LogNormalize", subtitle = "Descriptive only; no P values (amendment 01 section 5)"))
print(ggplot(ps, aes(sample, detection)) + geom_col() + facet_wrap(~population, ncol = 1) + theme_classic() +
        labs(y = "Fraction of cells with RAB7A UMI > 0", subtitle = "Descriptive only"))
dev.off()
write_json(list(script = "08_outcome_amendment01.R", run_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                review_record_sha256 = sha_file(rr_path), c10_sha256 = sha_file(c10_json),
                denominator = "colSums(RNA counts) == nCount_RNA; union of per-sample min.cells = 3 features",
                feature_universe = list(v1 = v1$feature_universe, v2 = v2$feature_universe),
                all_four_white_lower = as.list(all_lower), paired = pair,
                note = "Relative expression; no P values; v1 and v2 overlap and are not independent replications",
                session = capture.output(sessionInfo())),
           file.path(out, "E0_outcome.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
