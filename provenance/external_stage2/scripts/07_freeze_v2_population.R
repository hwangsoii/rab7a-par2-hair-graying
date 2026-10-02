# Amendment 01 (2026-10-01): freeze the v2 population (amendment section 4).
# Runs only if V2_D0 records conclusion == "provisionally_consistent" and every
# hash it binds (marker evidence, source object, report, script 06) is current.
# Membership is derived from the existing label; never reads RAB7A or F2RL1.
# Writes derived/stage2_v2_population_cells.txt (private), outputs/C10 (no
# self-checksum) and the UNAPPROVED review-record template. Run from the project root.
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "stage2_amend_common.R"))
plan_sha <- check_plan_frozen()
amend_sha <- check_amendment_frozen()
a4_sha <- check_a4_pass()
template <- file.path(stage2_dir, "REVIEW_RECORD_AMENDMENT01_TEMPLATE.json")
if (file.exists(c10_json) || file.exists(v2_cells_txt) || file.exists(template))
  stop("v2 already frozen; refusing to overwrite")

# Structured, hash-bound support decision
d0 <- fromJSON(d0_json, simplifyVector = FALSE)
if (!identical(d0$conclusion, "provisionally_consistent"))
  stop("V2_D0 conclusion is '", d0$conclusion, "': stop for review, no v2 freeze")
b <- d0$bound_sha256
bound <- c(list(list(name = "integrated_object", path = integrated_rds, sha256 = b$integrated_object),
                list(name = "C2_report", path = file.path(stage2_dir, "CHECKPOINT_C2_AMENDMENT01_MARKER_REVIEW.md"), sha256 = b$report_CHECKPOINT_C2),
                list(name = "script06", path = file.path(stage2_dir, "scripts", "06_marker_review_v2.R"), sha256 = b$script06),
                list(name = "amendment", path = amend_file, sha256 = b$amendment)),
           lapply(names(b$marker_evidence), function(f)
             list(name = f, path = file.path(out_dir, f), sha256 = b$marker_evidence[[f]])))
vb <- validate_dependencies(bound)
if (!vb$ok) { print(vb$failures); stop("V2_D0 binding broken") }

# v1 record still intact
c9 <- fromJSON(c9_json)
v1_cells <- readLines(v1_cells_txt)
stopifnot(identical(digest(paste(v1_cells, collapse = "\n"), algo = "sha256", serialize = FALSE), c9$population_cells_sha256),
          identical(sha_file(kc_rds), c9$kc_object_sha256))

# Derive v2 membership from the existing label
obj_sha <- sha_file(integrated_rds)
obj <- readRDS(integrated_rds)
md <- obj@meta.data; rm(obj); invisible(gc())
cand <- md$lineage == v2_label
cells <- sort(rownames(md)[cand])
per_sample <- table(factor(md$sample[cand], levels = samples$sample))
clusters <- sort(as.integer(as.character(unique(md$seurat_clusters[cand]))))
stopifnot(all(per_sample > 0), sum(per_sample) == length(cells), all(cells %in% v1_cells))
writeLines(cells, v2_cells_txt)

root <- normalizePath(getwd())
rel <- function(p) { p <- normalizePath(p); if (startsWith(p, paste0(root, "/"))) substring(p, nchar(root) + 2) else p }
dep <- function(name, path) list(name = name, path = rel(path), sha256 = sha_file(path))
sd <- file.path(stage2_dir, "scripts")
deps <- list(
  dep("plan", plan_file), dep("B0", freeze_txt), dep("C9", c9_json),
  dep("amendment", amend_file), dep("B1", amend_freeze_txt),
  dep("v1_cells", v1_cells_txt), dep("v2_cells", v2_cells_txt),
  dep("kc_object", kc_rds), dep("integrated_object", integrated_rds),
  dep("A4", a4_json), dep("V2_M0", file.path(out_dir, "V2_M0_marker_review.json")), dep("V2_D0", d0_json),
  dep("C2_report", file.path(stage2_dir, "CHECKPOINT_C2_AMENDMENT01_MARKER_REVIEW.md")),
  dep("script05", file.path(sd, "05_enforce_checkpoint_A.R")), dep("script06", file.path(sd, "06_marker_review_v2.R")),
  dep("script07", file.path(sd, "07_freeze_v2_population.R")), dep("script08", file.path(sd, "08_outcome_amendment01.R")),
  dep("stage2_common", file.path(sd, "stage2_common.R")), dep("stage2_validate", file.path(sd, "stage2_validate.R")),
  dep("stage2_amend_common", file.path(sd, "stage2_amend_common.R")))
stopifnot(identical(deps[[9]]$sha256, obj_sha))

write_json(list(record = "C10_v2_population_freeze", script = "07_freeze_v2_population.R",
                frozen_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                definitions = list(
                  v1 = "C9: all keratinocytes (rule_non_discriminating_all_keratinocytes); expression from kc object counts/data",
                  v2 = "all cells of stage2_integrated_annotated.rds with existing lineage == 'KC_cortex_cuticle'; expression from integrated object counts/data"),
                v2_derived_clusters_integrated_res0.5 = clusters, v2_n_cells = length(cells),
                v2_per_sample = as.list(per_sample), v2_subset_of_v1 = TRUE, v1_n_cells = length(v1_cells),
                annotation_decision = d0$conclusion,
                provisional_name = "hair-shaft differentiating keratinocytes with cortex/cuticle-marker enrichment",
                dependencies = deps,
                checksum_note = "C10 contains no checksum of itself; the independent review record stores C10's SHA-256.",
                session = capture.output(sessionInfo())),
           c10_json, auto_unbox = TRUE, pretty = TRUE, digits = NA)

# Unapproved review-record template (reviewers complete it; never pre-approved)
write_json(list(record = "independent review of Amendment 01 (v1 + v2 population freeze and runner 08)",
                approved = FALSE, reviewers = list(), review_date = NULL,
                c10_sha256 = sha_file(c10_json),
                approved_versions = setNames(lapply(deps, `[[`, "sha256"), vapply(deps, `[[`, "", "name")),
                instructions = paste("Reviewers verify each listed file and its SHA-256, then set approved = true,",
                                     "list at least two reviewer identities and an ISO review_date, and save as a NEW",
                                     "dated file. The template itself must remain approved = false.")),
           template, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat("v2 frozen:", length(cells), "cells; clusters", paste(clusters, collapse = ","), "\n")
print(per_sample)
