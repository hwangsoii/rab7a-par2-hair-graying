# Amendment 01 (2026-10-01): side-effect-free dependency validator.
# Reads files and computes SHA-256 only: no readRDS, no gene access, no
# directory creation. Sourced by stage2_amend_common.R and by the tests.
suppressPackageStartupMessages(library(digest))

sha_file <- function(path) digest(path, algo = "sha256", file = TRUE)

# entries: named list or data.frame of (name, path, sha256); paths are
# resolved against base_dir. Returns list(ok, failures).
validate_dependencies <- function(entries, base_dir = ".") {
  if (is.data.frame(entries)) entries <- split(entries, seq_len(nrow(entries)))
  fail <- list()
  add <- function(name, expected, observed, reason)
    fail[[length(fail) + 1]] <<- data.frame(name = name, expected = expected,
                                            observed = observed, reason = reason)
  for (i in seq_along(entries)) {
    e <- entries[[i]]
    name <- if (!is.null(e$name)) as.character(e$name) else names(entries)[i]
    if (is.null(e$path) || is.null(e$sha256) || is.na(e$path) || is.na(e$sha256) ||
        !nzchar(e$path) || !nzchar(e$sha256)) {
      add(name, NA_character_, NA_character_, "missing_field"); next
    }
    p <- if (startsWith(e$path, "/")) e$path else file.path(base_dir, e$path)
    if (!file.exists(p)) { add(name, e$sha256, NA_character_, "missing_file"); next }
    obs <- sha_file(p)
    if (!identical(obs, as.character(e$sha256))) add(name, e$sha256, obs, "hash_mismatch")
  }
  failures <- if (length(fail)) do.call(rbind, fail) else
    data.frame(name = character(), expected = character(), observed = character(), reason = character())
  list(ok = nrow(failures) == 0, failures = failures)
}
