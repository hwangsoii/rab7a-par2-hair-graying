# Amendment 01 (2026-10-01): fail-closed tests of the dependency validator and
# of runner 08's production approval gate. Validator tests use isolated
# test-only fixtures in <fixture_dir>; no approval record is used or created.
# Gate tests call 08 with --gate-only and NO genuine approval; each must stop
# with the expected reason before any object load or outcome-directory creation.
# The successful production-approval path is NOT tested (no genuine approval).
# Usage (from the project root): Rscript test_amendment01_gate.R <fixture_dir> <log_file>
args <- commandArgs(trailingOnly = TRUE)
fx <- args[1]; log_file <- args[2]
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(here, "..", "stage2_validate.R"))
stage2 <- normalizePath(file.path(here, "..", ".."))
res <- data.frame(test = character(), expected = character(), observed = character(), pass = logical())
rec <- function(test, expected, observed)
  res[nrow(res) + 1, ] <<- list(test, expected, paste(observed, collapse = "|"), identical(expected, paste(observed, collapse = "|")))

# ---- validator: isolated fixtures ------------------------------------------
dir.create(fx, recursive = TRUE, showWarnings = FALSE)
writeLines("fixture A", file.path(fx, "a.txt")); writeLines("fixture B", file.path(fx, "b.txt"))
man <- list(list(name = "a", path = "a.txt", sha256 = sha_file(file.path(fx, "a.txt"))),
            list(name = "b", path = "b.txt", sha256 = sha_file(file.path(fx, "b.txt"))))
v <- validate_dependencies(man, fx)
rec("validator_unchanged_ok", "TRUE", v$ok)
writeLines("fixture B.", file.path(fx, "b.txt"))  # one-byte change
v <- validate_dependencies(man, fx)
rec("validator_tampered_reason", "b:hash_mismatch", paste0(v$failures$name, ":", v$failures$reason))
file.remove(file.path(fx, "a.txt"))
v <- validate_dependencies(man[1], fx)
rec("validator_deleted_reason", "a:missing_file", paste0(v$failures$name, ":", v$failures$reason))
v <- validate_dependencies(list(list(name = "c", path = "c.txt")), fx)
rec("validator_missing_field_reason", "c:missing_field", paste0(v$failures$name, ":", v$failures$reason))

# ---- production gate: no approval / approved:false -------------------------
# Relative path: Rscript mangles spaces in an absolute --file path ("~+~")
runner <- file.path(here, "..", "08_outcome_amendment01.R")
outcome_dir <- file.path(stage2, "outputs", "E_outcome_amendment01")
gate <- function(extra) {
  o <- suppressWarnings(system2("Rscript", c(shQuote(runner), extra, "--gate-only"), stdout = TRUE, stderr = TRUE))
  list(status = if (is.null(attr(o, "status"))) 0L else attr(o, "status"), out = o)
}
pre_exists <- dir.exists(outcome_dir)
g <- gate(character())
rec("gate_no_record_reason", "GATE_FAIL: review_record_missing", trimws(grep("^GATE_FAIL", g$out, value = TRUE)))
rec("gate_no_record_exit", "2", g$status)
rec("gate_no_record_stages", "approval_record", trimws(sub("^GATE_STAGE: ", "", grep("^GATE_STAGE", g$out, value = TRUE))))
tmpl <- file.path(here, "..", "..", "REVIEW_RECORD_AMENDMENT01_TEMPLATE.json")
g <- gate(shQuote(paste0("--review-record=", tmpl)))
rec("gate_template_reason", "GATE_FAIL: review_not_approved", trimws(grep("^GATE_FAIL", g$out, value = TRUE)))
rec("gate_template_exit", "2", g$status)
rec("gate_template_stages", "approval_record", trimws(sub("^GATE_STAGE: ", "", grep("^GATE_STAGE", g$out, value = TRUE))))
rec("gate_never_ok", "0", length(grep("GATE_OK", g$out)))
rec("outcome_dir_not_created", "FALSE", pre_exists || dir.exists(outcome_dir))

out <- c(sprintf("Amendment 01 gate tests, %s", format(Sys.time(), tz = "UTC", usetz = TRUE)),
         sprintf("runner sha256: %s", sha_file(runner)),
         sprintf("validator sha256: %s", sha_file(file.path(here, "..", "stage2_validate.R"))),
         "Fixtures: isolated test-only files outside the project (no approval record used or created).",
         "Successful production-approval path: UNTESTED (no genuine approval exists).",
         capture.output(print(res, right = FALSE, row.names = FALSE)),
         sprintf("RESULT: %d/%d passed", sum(res$pass), nrow(res)))
writeLines(out, log_file); cat(out, sep = "\n")
if (!all(res$pass)) quit(status = 1)
