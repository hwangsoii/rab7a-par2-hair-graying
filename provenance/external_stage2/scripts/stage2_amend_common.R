# Amendment 01 (2026-10-01): shared paths, amendment-freeze check and marker
# lists. Sources the unchanged stage2_common.R and stage2_validate.R.
.amend_script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(.amend_script_dir, "stage2_common.R"))
source(file.path(.amend_script_dir, "stage2_validate.R"))

amend_file  <- file.path(stage2_dir, "AMENDMENT_01_POPULATION_v2_2026-10-01.md")
amend_freeze_txt <- file.path(out_dir, "B1_amendment01_freeze.txt")
a4_json     <- file.path(out_dir, "A4_checkpointA_enforced.json")
integrated_rds <- file.path(derived, "stage2_integrated_annotated.rds")
kc_rds      <- file.path(derived, "stage2_keratinocytes_subclustered.rds")
v1_cells_txt <- file.path(derived, "stage2_population_cells.txt")
v2_cells_txt <- file.path(derived, "stage2_v2_population_cells.txt")
c9_json     <- file.path(out_dir, "C9_population_freeze.json")
c10_json    <- file.path(out_dir, "C10_v2_population_freeze.json")
d0_json     <- file.path(out_dir, "V2_D0_annotation_decision.json")
v2_label    <- "KC_cortex_cuticle"

# Refuse to run unless the amendment's SHA-256 equals the one recorded in B1.
check_amendment_frozen <- function() {
  if (!file.exists(amend_freeze_txt)) stop("Amendment not frozen: ", amend_freeze_txt, " missing")
  rec <- sub("^amendment_sha256: ", "", grep("^amendment_sha256: ", readLines(amend_freeze_txt), value = TRUE))
  now <- sha_file(amend_file)
  if (!identical(rec, now)) stop("Amendment changed since freeze (", rec, " vs ", now, ")")
  cat("Amendment freeze verified:", now, "\n")
  invisible(now)
}

# Require a current, hash-bound A4 PASS (checkpoint-A enforcement).
check_a4_pass <- function() {
  if (!file.exists(a4_json)) stop("A4 missing: run 05_enforce_checkpoint_A.R first")
  a4 <- fromJSON(a4_json)
  if (!identical(a4$status, "PASS")) stop("A4 status is not PASS")
  bound <- list(A1 = file.path(out_dir, "A1_checksums.tsv"), A3 = file.path(out_dir, "A3_verification.json"),
                md5sum_txt = file.path(zen_dir, "md5sum.txt"),
                script05 = file.path(stage2_dir, "scripts", "05_enforce_checkpoint_A.R"))
  for (k in names(bound)) if (!identical(sha_file(bound[[k]]), a4$input_sha256[[k]]))
    stop("A4 binding broken for ", k)
  cat("A4 PASS verified\n")
  invisible(sha_file(a4_json))
}

# Annotation-review markers (amendment section 3). Outcome genes excluded.
review_markers <- list(
  cortex_cuticle = c("KRT31", "KRT35", "KRT85", "KRT32", "KRT33A"),
  basal_ORS = c("KRT5", "KRT14", "KRT15", "KRT17", "KRT6A", "KRT16", "COL17A1"),
  IRS = c("KRT25", "KRT27", "KRT71", "KRT73", "TCHH", "GATA3"),
  matrix_progenitor = c("MSX2", "LEF1", "HOXC13", "MKI67"),
  KRT1_KRT10 = c("KRT1", "KRT10"))
stopifnot(!any(outcome_genes %in% unlist(review_markers)))
