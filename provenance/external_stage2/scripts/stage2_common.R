# Shared paths and the plan-freeze gate for the Stage 2 scripts (2026-10-01).
suppressPackageStartupMessages({ library(digest); library(jsonlite) })

stage2_dir <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), ".."))
out_dir    <- file.path(stage2_dir, "outputs")
plan_file  <- file.path(stage2_dir, "ANALYSIS_PLAN_2026-10-01.md")
freeze_txt <- file.path(out_dir, "B0_plan_freeze.txt")
zen_dir    <- "OneDrive_2025-09-11/Data/cell discovery/onfroy_zenodo_15103193/download"
derived    <- "OneDrive_2025-09-11/Data/cell discovery/onfroy_zenodo_15103193/derived"
dir.create(derived, showWarnings = FALSE)

samples <- data.frame(
  sample    = c("F18", "F59", "F31B", "F31W", "F62B", "F62W"),
  fastq     = c("ryg035", "ryg029", "black-1", "white-2", "ryg047", "ryg048"),
  phenotype = c("Black", "Black", "Black", "White", "Black", "White"),
  donor     = c("F18", "F59", "F31", "F31", "F62", "F62"),
  stringsAsFactors = FALSE)

# Refuse to run unless the plan's SHA-256 equals the one recorded at freeze.
check_plan_frozen <- function() {
  if (!file.exists(freeze_txt)) stop("Plan not frozen: ", freeze_txt, " missing")
  rec <- sub("^plan_sha256: ", "", grep("^plan_sha256: ", readLines(freeze_txt), value = TRUE))
  now <- digest(plan_file, algo = "sha256", file = TRUE)
  if (!identical(rec, now)) stop("Plan changed since freeze (", rec, " vs ", now, ")")
  cat("Plan freeze verified:", now, "\n")
  invisible(now)
}

# Genes the pre-RAB7A scripts (02, 03) must never extract or report.
outcome_genes <- c("RAB7A", "F2RL1")
