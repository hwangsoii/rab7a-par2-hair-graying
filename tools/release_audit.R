#!/usr/bin/env Rscript
# =============================================================================
# release_audit.R — safety gate for a public code release tree
#
# Read-only. Never edits, deletes, stages, commits or uploads anything.
#
# Usage:
#   Rscript tools/release_audit.R [directory]
#
# Default directory is the repository root inferred from this script's location.
# Exit status 0 = pass, 1 = at least one FAIL.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

.script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grep("^--file=", a)]
  if (length(f) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", f[1])), mustWork = FALSE))
  }
  normalizePath(".", mustWork = FALSE)
}

root <- if (length(args) >= 1) args[1] else file.path(.script_dir(), "..")
root <- normalizePath(root, mustWork = TRUE)

cat("=============================================================\n")
cat("Release audit\n")
cat("Target:", root, "\n")
cat("=============================================================\n\n")

results <- list()
record <- function(check, status, detail = "") {
  results[[length(results) + 1]] <<- data.frame(
    check = check, status = status, detail = detail, stringsAsFactors = FALSE)
  cat(sprintf("[%-4s] %s%s\n", status, check,
              if (nzchar(detail)) paste0("\n         ", detail) else ""))
}

# --- File inventory ----------------------------------------------------------
all_files <- list.files(root, recursive = TRUE, all.files = TRUE,
                        full.names = TRUE, no.. = TRUE)
all_files <- all_files[!grepl("(^|/)\\.git/", all_files)]
all_files <- all_files[!dir.exists(all_files)]
rel <- sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", root), "/?"),
           "", all_files)

is_git_repo <- dir.exists(file.path(root, ".git"))
tracked <- NULL
if (is_git_repo) {
  tracked <- tryCatch(
    system2("git", c("-C", shQuote(root), "ls-files", "--cached", "--others",
                     "--exclude-standard"), stdout = TRUE),
    error = function(e) NULL)
}
# In a Git repo, audit tracked files plus non-ignored untracked files (everything
# currently eligible for a future commit); otherwise audit every file present in
# the candidate directory.
audit_rel <- if (!is.null(tracked) && length(tracked) > 0) tracked else rel
audit_abs <- file.path(root, audit_rel)
audit_abs <- audit_abs[file.exists(audit_abs)]
audit_rel <- sub(paste0("^", gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", root), "/?"),
                 "", audit_abs)

cat("Files audited:", length(audit_rel),
    if (!is.null(tracked) && length(tracked) > 0) {
      "(git-tracked + non-ignored untracked)"
    } else {
      "(all files)"
    },
    "\n\n")

# --- 1. R syntax -------------------------------------------------------------
r_files <- audit_abs[grepl("\\.[Rr]$", audit_abs)]
bad_parse <- character(0)
for (f in r_files) {
  ok <- tryCatch({ parse(f); TRUE }, error = function(e) {
    bad_parse <<- c(bad_parse, paste0(basename(f), ": ", conditionMessage(e)))
    FALSE })
}
if (length(bad_parse) == 0) {
  record("R syntax parses", "PASS", sprintf("%d R file(s)", length(r_files)))
} else {
  record("R syntax parses", "FAIL", paste(bad_parse, collapse = "\n         "))
}

# --- 2. Forbidden file extensions -------------------------------------------
forbidden_ext <- "\\.(fastq|fq|bam|sam|cram|h5ad|h5|h5seurat|loom|mtx|rds|rdata|rda|gz|zip|pptx|docx|ai|psd|tiff?|png|jpe?g|pdf|xlsx?|csv)$"
hits_ext <- audit_rel[grepl(forbidden_ext, audit_rel, ignore.case = TRUE)]
allowed_tabular <- c(
  "data/sample_manifest.tsv",
  "data/external/sample_manifest.tsv",
  "data/external/rab7a_cells.tsv",
  "data/external/sample_summaries.tsv",
  "data/external/paired_descriptive.tsv",
  "data/external/nominal_cell_tests.tsv",
  "data/external/cluster_selection.tsv",
  "data/external/v2_cluster_by_sample.tsv",
  "data/external/markers_by_cluster.tsv",
  "data/external/markers_by_lineage.tsv",
  "provenance/external_stage2/outputs/E_outcome_amendment01/E1_per_sample.tsv",
  "provenance/external_stage2/outputs/E_outcome_amendment01/E2_paired_donors.tsv"
)
hits_ext <- union(hits_ext, setdiff(
  audit_rel[grepl("\\.tsv$", audit_rel)], allowed_tabular))
if (length(hits_ext) == 0) {
  record("No data/binary/figure files", "PASS")
} else {
  record("No data/binary/figure files", "FAIL",
         paste(hits_ext, collapse = "\n         "))
}

# --- 3. Forbidden directories/documents -------------------------------------
forbidden_paths <- c("RAB7A-PAR2/", "OneDrive", "GEO_submission/",
                     "processed_data_release/", "release_staging/",
                     "config.local.R", "JID_release_audit", "Rab7A_JID",
                     "白", "논문")
hits_path <- audit_rel[Reduce(`|`, lapply(forbidden_paths, function(p)
  grepl(p, audit_rel, fixed = TRUE)))]
if (length(hits_path) == 0) {
  record("No manuscript/internal/deposition paths", "PASS")
} else {
  record("No manuscript/internal/deposition paths", "FAIL",
         paste(hits_path, collapse = "\n         "))
}

# --- Text scanning -----------------------------------------------------------
text_files <- audit_abs[grepl("\\.(R|r|py|json|md|txt|tsv|gitignore|yml|yaml)$|(^|/)\\.gitignore$",
                              audit_abs)]
read_text <- function(f) tryCatch(readLines(f, warn = FALSE),
                                  error = function(e) character(0))
scan_pattern <- function(pattern, label, ignore_case = TRUE, exclude_files = character(0)) {
  hits <- character(0)
  for (f in text_files) {
    if (basename(f) %in% exclude_files) next
    lines <- read_text(f)
    idx <- grep(pattern, lines, ignore.case = ignore_case, perl = TRUE)
    if (length(idx) > 0) {
      hits <- c(hits, sprintf("%s:%d: %s", sub(paste0("^", root, "/?"), "", f),
                              idx, trimws(substr(lines[idx], 1, 100))))
    }
  }
  hits
}

# --- 4. Absolute personal paths ---------------------------------------------
# Fragments again, so this file does not match its own patterns.
path_tokens <- c(
  paste0("/Us", "ers/[A-Za-z0-9._-]+"),
  paste0("/ho", "me/[A-Za-z0-9._-]+"),
  paste0("[A-Z]:\\\\\\\\Us", "ers"),
  paste0("Cloud", "Storage"),
  paste0("Google", "Drive-"),
  paste0("My ", "Drive")
)
hits_abs <- scan_pattern(paste0("(", paste(path_tokens, collapse = "|"), ")"),
                         "abs paths")
if (length(hits_abs) == 0) {
  record("No absolute or cloud-account paths", "PASS")
} else {
  record("No absolute or cloud-account paths", "FAIL",
         paste(utils::head(hits_abs, 10), collapse = "\n         "))
}

# --- 5. E-mail addresses and credentials ------------------------------------
# The e-mail pattern requires a real top-level domain so that R slot access
# (obj@meta.data) is not mistaken for an address.
tld <- "(com|org|net|edu|gov|mil|int|io|co|ac|kr|uk|de|jp|cn|fr|it|es|nl|se|ch|info|mail)"
email_re <- paste0("[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\\.", tld, "\\b")
hits_cred <- scan_pattern(paste0(email_re,
                                 "|password\\s*[:=]|secret\\s*[:=]|api[_-]?key\\s*[:=]",
                                 "|BEGIN [A-Z ]*PRIVATE KEY|ghp_[A-Za-z0-9]{20,}"),
                          "credentials")
if (length(hits_cred) == 0) {
  record("No e-mail addresses or credentials", "PASS")
} else {
  record("No e-mail addresses or credentials", "FAIL",
         paste(utils::head(hits_cred, 10), collapse = "\n         "))
}

# --- 6. Prohibited participant metadata -------------------------------------
# Column headers and value patterns that must not appear in released metadata.
tsv_files <- audit_abs[grepl("\\.tsv$", audit_abs)]
bad_cols <- character(0)
prohibited <- c("age", "sex", "gender", "birth", "dob", "collection_year",
                "collection_date", "date", "initials", "patient", "mrn",
                "hospital", "institution", "ethnic", "race")
for (f in tsv_files) {
  hdr <- tryCatch(strsplit(readLines(f, n = 1), "\t")[[1]],
                  error = function(e) character(0))
  hit <- hdr[tolower(hdr) %in% prohibited]
  if (length(hit) > 0) {
    bad_cols <- c(bad_cols, sprintf("%s: %s",
                                    sub(paste0("^", root, "/?"), "", f),
                                    paste(hit, collapse = ", ")))
  }
}
if (length(bad_cols) == 0) {
  record("No prohibited metadata columns", "PASS",
         sprintf("%d tsv file(s) checked", length(tsv_files)))
} else {
  record("No prohibited metadata columns", "FAIL",
         paste(bad_cols, collapse = "\n         "))
}

# --- 7. IRB / institution / donor identifiers in text ------------------------
# Tokens are assembled from fragments so that this audit script does not itself
# contain the identifiers it searches for.
sensitive_tokens <- c(
  paste0("IRB", " No"),
  paste0("4-", "2022-", "0778"),
  paste0("Sever", "ance"),
  paste0("Yon", "sei"),
  paste0("\\bO", "SH\\b"),
  paste0("\\bT", "GK\\b")
)
hits_irb <- scan_pattern(paste(sensitive_tokens, collapse = "|"),
                         "institutional identifiers")
if (length(hits_irb) == 0) {
  record("No institutional or donor identifiers", "PASS")
} else {
  record("No institutional or donor identifiers", "FAIL",
         paste(utils::head(hits_irb, 10), collapse = "\n         "))
}

# --- 8. Fabricated or stale accessions --------------------------------------
acc_tokens <- c(paste0("GS", "E[X0-9]{4,}"),
                paste0("SR", "[PRXS][X0-9]{4,}"),
                paste0("PRJ", "[EDN]A[X0-9]{3,}"))
hits_acc <- scan_pattern(paste(acc_tokens, collapse = "|"), "accessions")
if (length(hits_acc) == 0) {
  record("No claimed/placeholder accessions", "PASS")
} else {
  record("No claimed/placeholder accessions", "FAIL",
         paste(utils::head(hits_acc, 10), collapse = "\n         "))
}

# --- 9. No unresolved placeholders ------------------------------------------
# Assemble tokens from fragments so the audit definition is not counted as a
# placeholder occurrence.
placeholder_tokens <- c(
  paste0("\\[PROCESSED_DATA_", "REPOSITORY_LINK\\]"),
  paste0("\\[PUBLIC_CODE_", "REPOSITORY_LINK\\]"),
  paste0("<[A-Z_]+_", "TBD>"),
  paste0("<COPYRIGHT_HOLDER_", "NAME>")
)
hits_ph <- scan_pattern(paste(placeholder_tokens, collapse = "|"),
                        "placeholders")
if (length(hits_ph) == 0) {
  record("No unresolved placeholders", "PASS")
} else {
  record("No unresolved placeholders", "FAIL",
         paste(utils::head(hits_ph, 10), collapse = "\n         "))
}

# --- 10. Required documentation ---------------------------------------------
required_docs <- c("README.md", "LICENSE", "DATA_DICTIONARY.md",
                   "REPRODUCIBILITY.md", "PROCESSED_DATA_PACKAGE_README.md",
                   "supplementary_code/README.md")
missing_docs <- required_docs[!file.exists(file.path(root, required_docs))]
if (length(missing_docs) == 0) {
  record("Required documentation present", "PASS")
} else {
  record("Required documentation present", "FAIL",
         paste(missing_docs, collapse = ", "))
}

# --- 11. Broken internal documentation links --------------------------------
md_files <- audit_abs[grepl("\\.md$", audit_abs)]
broken <- character(0)
for (f in md_files) {
  lines <- read_text(f)
  links <- unlist(regmatches(lines, gregexpr("\\]\\(([^)]+)\\)", lines)))
  links <- gsub("^\\]\\(|\\)$", "", links)
  links <- links[!grepl("^(https?:|mailto:|#)", links)]
  for (l in links) {
    target <- file.path(dirname(f), sub("#.*$", "", l))
    if (!file.exists(target) && !dir.exists(target)) {
      broken <- c(broken, sprintf("%s -> %s",
                                  sub(paste0("^", root, "/?"), "", f), l))
    }
  }
}
if (length(broken) == 0) {
  record("Internal documentation links resolve", "PASS")
} else {
  record("Internal documentation links resolve", "FAIL",
         paste(utils::head(broken, 10), collapse = "\n         "))
}

# --- 12. Git presence in a release candidate --------------------------------
if (is_git_repo) {
  record("Git repository present", "INFO",
         "internal repository — clean release candidates must have no .git")
} else {
  record("No .git directory", "PASS", "clean release candidate")
}

# --- Summary -----------------------------------------------------------------
summary_df <- do.call(rbind, results)
n_fail <- sum(summary_df$status == "FAIL")
cat("\n-------------------------------------------------------------\n")
cat(sprintf("PASS: %d   FAIL: %d   INFO: %d\n",
            sum(summary_df$status == "PASS"), n_fail,
            sum(summary_df$status == "INFO")))
cat("-------------------------------------------------------------\n")
if (n_fail > 0) {
  cat("Release audit FAILED. Do not publish this tree.\n")
  quit(status = 1)
}
cat("Release audit passed for the checks implemented here.\n")
cat("Passing is necessary, not sufficient: corresponding-author review and\n")
cat("editorial approval are still required.\n")
quit(status = 0)
