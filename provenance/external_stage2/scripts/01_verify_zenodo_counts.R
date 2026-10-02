# Stage 2, checkpoint A (2026-10-01): verify the third-party count matrices
# (Zenodo 15103193, Onfroy 2025) before any analysis.
# Questions: checksums; genuine nonnegative integer UMI counts; raw-droplet vs
# Cell Ranger-filtered matrices; per-sample Cell Ranger records; consistency of
# the sample labels with the official NODE FASTQ mapping; and whether the local
# h5ad sample names correspond to these samples (barcode overlap, verification
# only -- the old objects are NOT analysis inputs).
# No expression value of any individual gene is extracted. Sparse only; one
# sample at a time. Writes aggregate JSON/TSV only (no barcodes).
# Run from the project root.
suppressPackageStartupMessages({ library(Matrix); library(jsonlite); library(digest) })

zen_dir <- "OneDrive_2025-09-11/Data/cell discovery/onfroy_zenodo_15103193/download"
local_meta <- "OneDrive_2025-09-11/Data/cell discovery/seurat_data/metadata.tsv"
manifest_xlsx <- "OneDrive_2025-09-11/Data/cell discovery/OEP00002321_Data_1760192653404.xlsx"
out_dir <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                     "..", "outputs")
samples <- c("F18", "F59", "F31B", "F31W", "F62B", "F62W")
res <- list(audit_date = "2026-10-01", script = "01_verify_zenodo_counts.R")

# ---- 1. checksums ----------------------------------------------------------
rec <- fromJSON(file.path(zen_dir, "zenodo_record_15103193.json"))
files <- rec$files
files$md5_expected <- sub("^md5:", "", files$checksum)
files$md5_observed <- vapply(file.path(zen_dir, files$key), function(f)
  if (file.exists(f)) digest(f, algo = "md5", file = TRUE) else NA_character_, "")
files$size_observed <- file.size(file.path(zen_dir, files$key))
files$md5_match <- files$md5_expected == files$md5_observed
# md5sum.txt shipped inside the record (second, author-side list)
md5txt <- read.table(file.path(zen_dir, "md5sum.txt"), stringsAsFactors = FALSE,
                     col.names = c("md5", "path"))
md5txt$key <- basename(md5txt$path)
res$checksums <- list(
  n_files = nrow(files), n_md5_match_zenodo_api = sum(files$md5_match, na.rm = TRUE),
  all_sizes_match = all(files$size_observed == files$size),
  md5sum_txt_rows = nrow(md5txt), md5sum_txt_paths = md5txt$path)
write.table(files[, c("key", "size", "md5_expected", "md5_observed", "md5_match")],
            file.path(out_dir, "A1_checksums.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(all(files$md5_match))

# ---- 2. per-sample matrix checks (sequential, sparse) ------------------------
mtx_header <- function(f) { con <- gzfile(f); on.exit(close(con)); readLines(con, n = 1) }
read_one <- function(s, kind) {
  p <- function(x) file.path(zen_dir, sprintf("%s_%s_%s", s, kind, x))
  m <- as(readMM(gzfile(p("matrix.mtx.gz"))), "CsparseMatrix")
  feats <- read.delim(gzfile(p("features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  bcs <- readLines(gzfile(p("barcodes.tsv.gz")))
  stopifnot(nrow(m) == nrow(feats), ncol(m) == length(bcs))
  list(m = m, feats = feats, bcs = bcs, header = mtx_header(p("matrix.mtx.gz")))
}
web_field <- function(s, field) {
  t <- readLines(file.path(zen_dir, sprintf("%s_web_summary.html", s)), warn = FALSE)
  t <- paste(t, collapse = "\n")
  m <- regmatches(t, regexpr(sprintf('\\["%s", "[^"]*"\\]', field), t))
  if (length(m)) sub('.*", "([^"]*)"\\]$', "\\1", m) else NA_character_
}
ref_feats <- NULL
per_sample <- list()
for (s in samples) {
  cat("sample", s, "\n")
  raw <- read_one(s, "raw"); flt <- read_one(s, "filtered")
  if (is.null(ref_feats)) ref_feats <- raw$feats
  in_raw <- match(flt$bcs, raw$bcs)
  sub_raw <- raw$m[, in_raw[!is.na(in_raw)], drop = FALSE]
  same_vals <- !anyNA(in_raw) && identical(dim(sub_raw), dim(flt$m)) && {
    d <- drop0(sub_raw - flt$m); length(d@x) == 0 }
  cs_flt <- colSums(flt$m); cs_raw <- colSums(raw$m)
  nf_flt <- diff(flt$m@p)
  bc_len <- unique(nchar(sub("-1$", "", flt$bcs)))
  per_sample[[s]] <- list(
    raw = list(dim = dim(raw$m), nnz = length(raw$m@x), mtx_header = raw$header,
               all_nonneg_integer_FULL = all(raw$m@x >= 0 & raw$m@x == round(raw$m@x)),
               min_nonzero = min(raw$m@x), max = max(raw$m@x),
               n_barcodes_with_any_umi = sum(cs_raw > 0)),
    filtered = list(dim = dim(flt$m), nnz = length(flt$m@x), mtx_header = flt$header,
                    all_nonneg_integer_FULL = all(flt$m@x >= 0 & flt$m@x == round(flt$m@x)),
                    min_nonzero = min(flt$m@x), max = max(flt$m@x),
                    umi_per_barcode = as.list(setNames(quantile(cs_flt, c(0, .5, 1)), c("min", "median", "max"))),
                    genes_per_barcode_median = median(nf_flt),
                    barcode_suffix_all_dash1 = all(grepl("-1$", flt$bcs)), barcode_core_length = bc_len),
    filtered_subset_of_raw = !anyNA(in_raw),
    filtered_values_identical_to_raw_columns_FULL = isTRUE(same_vals),
    filtered_are_top_raw_by_umi = min(cs_flt) >= stats::quantile(cs_raw[-in_raw], 1),  # informative only
    features_identical_raw_vs_filtered = identical(raw$feats, flt$feats),
    features_identical_to_F18 = identical(raw$feats, ref_feats),
    feature_types = as.list(table(raw$feats$V3)),
    web_summary = list(sample_id = web_field(s, "Sample ID"), chemistry = web_field(s, "Chemistry"),
                       transcriptome = web_field(s, "Transcriptome"),
                       pipeline_version = web_field(s, "Pipeline Version"),
                       number_of_reads = web_field(s, "Number of Reads"),
                       estimated_cells = web_field(s, "Estimated Number of Cells"),
                       median_umi = web_field(s, "Median UMI Counts per Cell"),
                       median_genes = web_field(s, "Median Genes per Cell")),
    filtered_ncells_equals_web_estimate =
      ncol(flt$m) == as.integer(gsub(",", "", web_field(s, "Estimated Number of Cells"))),
    filtered_median_umi_equals_web =
      median(cs_flt) == as.numeric(gsub(",", "", web_field(s, "Median UMI Counts per Cell"))))
  # keep only barcode cores for the identity check, in memory
  assign(paste0("bc_", s), sub("-1$", "", flt$bcs))
  assign(paste0("umi_", s), setNames(cs_flt, sub("-1$", "", flt$bcs)))
  assign(paste0("ngene_", s), setNames(nf_flt, sub("-1$", "", flt$bcs)))
  rm(raw, flt, sub_raw); invisible(gc())
}
res$per_sample <- per_sample
res$n_features <- nrow(ref_feats)
res$key_genes_present <- as.list(setNames(c("RAB7A", "KRT35", "KRT85", "KRT31", "F2RL1") %in% ref_feats$V2,
                                          c("RAB7A", "KRT35", "KRT85", "KRT31", "F2RL1")))
res$key_genes_unique <- all(table(ref_feats$V2)[c("RAB7A", "KRT35", "KRT85")] == 1)

# ---- 3. FASTQ-side evidence from the official NODE manifest -----------------
# Read via python-free route: the manifest is xlsx; extracted by 01b helper
# (sizes/run names only) into a TSV beforehand.
man <- read.delim(file.path(out_dir, "A0_node_manifest_by_prefix.tsv"), stringsAsFactors = FALSE)
res$node_manifest_by_prefix <- man

# ---- 4. local h5ad sample names vs Zenodo samples (verification only) -------
meta <- read.delim(local_meta, stringsAsFactors = FALSE, check.names = FALSE)
colnames(meta)[1] <- "cell"
meta$core <- sub("-1_[0-9]+$", "", meta$cell)
ov <- do.call(rbind, lapply(unique(meta$sample), function(ls) {
  m <- meta[meta$sample == ls, ]
  do.call(rbind, lapply(samples, function(s) {
    bc <- get(paste0("bc_", s)); hit <- m$core %in% bc
    umi <- get(paste0("umi_", s))[m$core[hit]]
    ng <- get(paste0("ngene_", s))[m$core[hit]]
    data.frame(local_sample = ls, zenodo_sample = s, local_cells = nrow(m),
               zenodo_cells = length(bc), overlap = sum(hit),
               frac_local_in_zenodo = round(mean(hit), 4),
               total_counts_equal = if (any(hit)) sum(umi == m$total_counts[hit]) else 0,
               n_genes_equal = if (any(hit)) sum(ng == m$n_genes_by_counts[hit]) else 0,
               median_ratio_total_counts = if (any(hit)) round(median(m$total_counts[hit] / umi), 4) else NA)
  }))
}))
write.table(ov, file.path(out_dir, "A2_local_vs_zenodo_barcode_overlap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
best <- do.call(rbind, lapply(split(ov, ov$local_sample), function(d) d[which.max(d$overlap), ]))
res$local_name_link <- best

dir.create(out_dir, showWarnings = FALSE)
write_json(res, file.path(out_dir, "A3_verification.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
print(best)
cat(toJSON(lapply(per_sample, function(x) x[c("filtered_subset_of_raw",
  "filtered_values_identical_to_raw_columns_FULL", "filtered_ncells_equals_web_estimate",
  "filtered_median_umi_equals_web")]), auto_unbox = TRUE, pretty = TRUE))
