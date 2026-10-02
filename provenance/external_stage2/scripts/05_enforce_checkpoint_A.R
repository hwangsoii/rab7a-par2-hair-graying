# Amendment 01 (2026-10-01): fail-closed enforcement of the RECORDED
# checkpoint-A results (A1, A3) from 01_verify_zenodo_counts.R. This wrapper
# does not repeat the full matrix checks; it asserts every required saved
# criterion, joins the author's md5sum.txt to the observed MD5s (42/42 hash
# equality), and re-hashes the 43 downloaded files to show they are unchanged.
# Writes outputs/A4_checkpointA_enforced.json. Run from the project root.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_amend_common.R"))

a1_tsv <- file.path(out_dir, "A1_checksums.tsv")
a3 <- file.path(out_dir, "A3_verification.json")
md5txt_file <- file.path(zen_dir, "md5sum.txt")
script01 <- file.path(stage2_dir, "scripts", "01_verify_zenodo_counts.R")
self <- file.path(stage2_dir, "scripts", "05_enforce_checkpoint_A.R")

checks <- list()
req <- function(name, value) {
  ok <- length(value) == 1 && !is.na(value) && isTRUE(as.logical(value))
  checks[[length(checks) + 1]] <<- data.frame(check = name, ok = ok)
  invisible(ok)
}
field <- function(x, ...) { for (k in c(...)) { if (is.null(x[[k]])) return(NULL); x <- x[[k]] }; x }

# Verifier version: script 01 must be the version recorded at plan freeze (B0)
b0 <- readLines(freeze_txt)
sha01_b0 <- sub("^\\s*([0-9a-f]{64})\\s+01_verify_zenodo_counts.R$", "\\1",
                grep("01_verify_zenodo_counts.R$", b0, value = TRUE))
req("verifier_01_sha_equals_B0", identical(sha_file(script01), sha01_b0))

A <- fromJSON(a3, simplifyVector = FALSE)
a1 <- read.delim(a1_tsv, stringsAsFactors = FALSE)
req("A3_checksums_n_files_43", identical(as.integer(field(A, "checksums", "n_files")), 43L))
req("A3_md5_match_zenodo_api_43", identical(as.integer(field(A, "checksums", "n_md5_match_zenodo_api")), 43L))
req("A3_all_sizes_match", field(A, "checksums", "all_sizes_match"))
req("A1_rows_43", nrow(a1) == 43)
req("A1_all_md5_match", all(a1$md5_match %in% TRUE) && !anyNA(a1$md5_observed))

for (s in samples$sample) {
  ps <- field(A, "per_sample", s)
  req(paste0(s, "_present"), !is.null(ps))
  for (kind in c("raw", "filtered")) {
    req(paste(s, kind, "all_nonneg_integer_FULL"), field(ps, kind, "all_nonneg_integer_FULL"))
    req(paste(s, kind, "header_coordinate_integer"),
        identical(field(ps, kind, "mtx_header"), "%%MatrixMarket matrix coordinate integer general"))
    req(paste(s, kind, "min_nonzero_ge_1"), isTRUE(field(ps, kind, "min_nonzero") >= 1))
    req(paste(s, kind, "n_features_27955"), identical(as.integer(field(ps, kind, "dim")[[1]]), 27955L))
  }
  req(paste(s, "filtered_subset_of_raw"), field(ps, "filtered_subset_of_raw"))
  req(paste(s, "filtered_values_identical_to_raw"), field(ps, "filtered_values_identical_to_raw_columns_FULL"))
  req(paste(s, "features_identical_raw_vs_filtered"), field(ps, "features_identical_raw_vs_filtered"))
  req(paste(s, "features_identical_to_F18"), field(ps, "features_identical_to_F18"))
  ft <- field(ps, "feature_types")
  req(paste(s, "feature_types_all_gene_expression"),
      identical(names(ft), "Gene Expression") && identical(as.integer(ft[[1]]), 27955L))
  req(paste(s, "web_sample_id_equals"), identical(field(ps, "web_summary", "sample_id"), s))
  req(paste(s, "web_ncells_equal"), field(ps, "filtered_ncells_equals_web_estimate"))
  req(paste(s, "web_median_umi_equal"), field(ps, "filtered_median_umi_equals_web"))
}
req("n_features_27955", identical(as.integer(field(A, "n_features")), 27955L))
for (g in c("RAB7A", "KRT35", "KRT85", "KRT31", "F2RL1")) req(paste0("gene_present_", g), field(A, "key_genes_present", g))
req("key_genes_unique", field(A, "key_genes_unique"))

# Author manifest: actual hash comparison, not a row count
m <- read.table(md5txt_file, stringsAsFactors = FALSE, col.names = c("md5", "path"))
m$key <- basename(m$path)
j <- merge(m, a1, by = "key", all.x = TRUE)
req("author_manifest_rows_42", nrow(m) == 42)
req("author_manifest_all_keys_in_A1", !anyNA(j$md5_observed))
req("author_manifest_42_of_42_hash_equal", sum(j$md5 == j$md5_observed, na.rm = TRUE) == 42)

# Current integrity: re-hash the 43 downloaded files against the Zenodo MD5s
now_md5 <- vapply(file.path(zen_dir, a1$key), function(f)
  if (file.exists(f)) digest(f, algo = "md5", file = TRUE) else NA_character_, "")
req("download_43_files_unchanged_md5", all(now_md5 == a1$md5_expected, na.rm = FALSE) && !anyNA(now_md5))
req("download_43_files_size", all(file.size(file.path(zen_dir, a1$key)) == a1$size))

checks <- do.call(rbind, checks)
status <- if (all(checks$ok)) "PASS" else "FAIL"
ws <- lapply(samples$sample, function(s) field(A, "per_sample", s, "web_summary")[c("pipeline_version", "transcriptome")])
names(ws) <- samples$sample
write_json(list(script = "05_enforce_checkpoint_A.R", run_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                status = status, n_checks = nrow(checks), n_failed = sum(!checks$ok),
                failed = checks$check[!checks$ok], checks = checks,
                scope = paste("Enforces results recorded by 01_verify_zenodo_counts.R (A1, A3);",
                              "does not independently repeat the full raw/filtered matrix checks.",
                              "Author-manifest comparison is a 42/42 hash join; downloaded files re-hashed (MD5)."),
                input_sha256 = list(A1 = sha_file(a1_tsv), A3 = sha_file(a3), md5sum_txt = sha_file(md5txt_file),
                                    script01 = sha_file(script01), script05 = sha_file(self)),
                known_discrepancies = list(
                  note = paste("Machine records report the pipeline/reference below; deposit prose says",
                               "Cell Ranger 3.0.1 / hg19. Recorded, not a failure; effect on quantification untested."),
                  web_summary = ws),
                reviewer_A_receipt = paste("Separate record (Codex checkpoint-C review): 43/43 checksum-and-size,",
                                           "42/42 author-manifest hashes, 36 barcode pairings, 12 matrix headers, XLSX mapping."),
                session = capture.output(sessionInfo())),
           a4_json, auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("A4 status:", status, "-", nrow(checks), "checks,", sum(!checks$ok), "failed\n")
if (status != "PASS") { print(checks[!checks$ok, ]); quit(status = 1) }
