# Downstream-only reproduction of the two frozen count-derived populations.
# Requires ggplot2, digest, jsonlite; no Seurat objects or network access.
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat("Rscript 01_make_Figure_S2B.R [--data-dir DIR] [--output-dir NEW_DIR] [--font-family FONT]\n")
  quit(save = "no", status = 0)
}
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE))
settings <- list(
  "data-dir" = file.path(script_dir, "..", "..", "data", "external"),
  "output-dir" = file.path(getwd(), "external_cohort_replot"), "font-family" = "Arial"
)
if (length(args) %% 2L != 0L) stop("Options require a value; see --help")
seen <- character()
for (i in (seq_len(length(args) %/% 2L) * 2L - 1L)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || !key %in% names(settings) || key %in% seen) {
    stop("Unknown or repeated option: ", args[i])
  }
  if (!nzchar(args[i + 1L])) stop("Empty option value")
  settings[[key]] <- args[i + 1L]
  seen <- c(seen, key)
}
for (package in c("ggplot2", "digest", "jsonlite")) {
  if (!requireNamespace(package, quietly = TRUE)) stop(package, " is required; nothing is installed automatically")
}
suppressPackageStartupMessages(library(ggplot2))
if (!capabilities("cairo")) stop("This reproduction requires R with Cairo PDF/PNG support")
data_dir <- normalizePath(settings[["data-dir"]], mustWork = TRUE)
output_dir <- path.expand(settings[["output-dir"]])
if (file.exists(output_dir)) stop("Output directory already exists; choose a new directory")
font_family <- settings[["font-family"]]
table_names <- c("sample_manifest", "rab7a_cells", "sample_summaries", "paired_descriptive",
                 "nominal_cell_tests", "cluster_selection", "v2_cluster_by_sample",
                 "markers_by_cluster", "markers_by_lineage")
input_paths <- file.path(data_dir, paste0(table_names, ".tsv"))
stopifnot(all(file.exists(input_paths)))
sha_file <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
export_qa <- jsonlite::fromJSON(file.path(data_dir, "export_qa.json"))
expected_paths <- file.path("data", "external", paste0(table_names, ".tsv"))
stopifnot(nrow(export_qa$files) == length(table_names),
          !anyDuplicated(export_qa$files$path), setequal(export_qa$files$path, expected_paths))
for (i in seq_along(input_paths)) {
  expected <- export_qa$files$sha256[match(expected_paths[i], export_qa$files$path)]
  if (sha_file(input_paths[i]) != expected) stop("Package source checksum mismatch: ", basename(input_paths[i]))
}
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
if (sha_file(script_path) != export_qa$replot_script$sha256) stop("Package replot script checksum mismatch")
input_hashes <- tools::md5sum(input_paths)
read_tsv <- function(name) {
  value <- read.delim(file.path(data_dir, paste0(name, ".tsv")), check.names = FALSE,
                      stringsAsFactors = FALSE,
                      colClasses = if (name == "rab7a_cells") c(rab7a_lognorm_hex = "character") else NA)
  stopifnot(!anyNA(value))
  value
}
samples <- read_tsv("sample_manifest")
cells <- read_tsv("rab7a_cells")
# Some R decimal readers shift a last binary digit; hex preserves source ties exactly.
stopifnot(all(grepl("^0x[0-9a-f]+(\\.[0-9a-f]+)?p[+-][0-9]+$", cells$rab7a_lognorm_hex)))
exact_expression <- as.numeric(cells$rab7a_lognorm_hex)
stopifnot(all(is.finite(exact_expression)),
          max(abs(cells$rab7a_lognorm - exact_expression)) < 1e-14)
cells$rab7a_lognorm <- exact_expression
summaries <- read_tsv("sample_summaries")
pair_reference <- read_tsv("paired_descriptive")
test_reference <- read_tsv("nominal_cell_tests")
selection <- read_tsv("cluster_selection")
cluster_counts <- read_tsv("v2_cluster_by_sample")
stopifnot(
  identical(names(cells), c("row_id", "population", "sample", "phenotype", "rab7a_umi", "rab7a_lognorm", "rab7a_lognorm_hex")),
  identical(cells$row_id, sprintf("E%06d", seq_len(nrow(cells)))),
  !anyDuplicated(cells$row_id), nrow(cells) == 18529L,
  identical(samples$sample, c("F18", "F59", "F31B", "F31W", "F62B", "F62W")),
  identical(samples$public_fastq_prefix, c("ryg035", "ryg029", "black-1", "white-2", "ryg047", "ryg048")),
  identical(samples$phenotype, c("black", "black", "black", "white", "black", "white")),
  identical(samples$donor_group, c("D01", "D02", "D03", "D03", "D04", "D04")),
  all(cells$sample %in% samples$sample),
  identical(cells$phenotype, samples$phenotype[match(cells$sample, samples$sample)]),
  all(is.finite(cells$rab7a_lognorm)), all(cells$rab7a_lognorm >= 0),
  all(is.finite(cells$rab7a_umi)), all(cells$rab7a_umi >= 0),
  all(cells$rab7a_umi == round(cells$rab7a_umi)),
  all((cells$rab7a_lognorm == 0) == (cells$rab7a_umi == 0)),
  nrow(summaries) == 12L, !anyDuplicated(summaries[c("population", "sample")]),
  nrow(test_reference) == 2L, nrow(pair_reference) == 4L,
  sum(selection$n_cells[selection$in_v1]) == 15214L,
  sum(selection$n_cells[selection$in_v2]) == 3315L,
  setequal(selection$cluster[selection$in_v2], c(3L, 6L, 9L, 14L)),
  nrow(cluster_counts) == 24L, sum(cluster_counts$n_cells) == 3315L
)
definitions <- list(
  list(id = "v2_KC_cortex_cuticle", title = "KRT35/KRT85-enriched\nkeratinocytes",
       filename = "External_RAB7A_differentiating_KC_cell_level", n = c(2684L, 631L),
       sample_n = c(179L, 820L, 547L, 247L, 1138L, 384L), w = 952884, p = 8.0794986333324072e-07),
  list(id = "v1_all_keratinocytes", title = "All keratinocytes",
       filename = "External_RAB7A_all_KC_cell_level", n = c(10382L, 4832L),
       sample_n = c(1436L, 2539L, 3012L, 2487L, 3395L, 2345L), w = 26497275.5, p = 9.4861419954486193e-09)
)
stopifnot(setequal(cells$population, vapply(definitions, `[[`, character(1), "id")),
          all(test_reference$alternative == "two.sided"), !any(test_reference$paired),
          !any(test_reference$exact), all(test_reference$correct),
          !any(test_reference$donor_adjusted), !any(test_reference$multiple_testing_adjusted))
sample_differences <- numeric()
prepared <- lapply(definitions, function(definition) {
  values <- cells[cells$population == definition$id, ]
  values$phenotype <- factor(values$phenotype, levels = c("black", "white"))
  stopifnot(identical(as.integer(table(values$phenotype)), definition$n))
  for (j in seq_len(nrow(samples))) {
    selected <- values$sample == samples$sample[j]
    reference <- summaries[summaries$population == definition$id & summaries$sample == samples$sample[j], ]
    stopifnot(nrow(reference) == 1L, reference$n_cells == definition$sample_n[j],
              sum(selected) == reference$n_cells,
              sum(values$rab7a_umi[selected]) == reference$rab7a_umi,
              sum(values$rab7a_umi[selected] == 0) == reference$n_zero,
              is.finite(reference$denominator_umi), reference$denominator_umi > 0,
              reference$phenotype == samples$phenotype[j], reference$donor_group == samples$donor_group[j])
    differences <- c(
      abs(mean(values$rab7a_lognorm[selected]) - reference$mean_lognorm),
      abs(mean(values$rab7a_umi[selected] > 0) - reference$detection),
      abs(log2(reference$rab7a_umi / reference$denominator_umi * 1e6 + 1) - reference$pseudobulk_log2cpm1)
    )
    stopifnot(max(differences) < 1e-12)
    sample_differences <<- c(sample_differences, differences)
  }
  result <- stats::wilcox.test(
    values$rab7a_lognorm[values$phenotype == "black"],
    values$rab7a_lognorm[values$phenotype == "white"],
    alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE
  )
  reference <- test_reference[test_reference$population == definition$id, ]
  stopifnot(nrow(reference) == 1L, unname(result$statistic) == definition$w,
            abs(result$p.value / definition$p - 1) < 1e-12,
            unname(result$statistic) == reference$statistic_w,
            abs(result$p.value / reference$p_nominal - 1) < 1e-12,
            reference$n_black == definition$n[1], reference$n_white == definition$n[2],
            sum(values$rab7a_lognorm == 0 & values$phenotype == "black") == reference$zeros_black,
            sum(values$rab7a_lognorm == 0 & values$phenotype == "white") == reference$zeros_white)
  list(definition = definition, values = values, result = result,
       p_label = paste0("Nominal P = ", format.pval(result$p.value, digits = 3, eps = 1e-300)))
})

paired <- do.call(rbind, lapply(definitions, function(definition) {
  do.call(rbind, lapply(c("D03", "D04"), function(donor_group) {
    pair <- summaries[summaries$population == definition$id & summaries$donor_group == donor_group, ]
    black <- pair[pair$phenotype == "black", ]
    white <- pair[pair$phenotype == "white", ]
    stopifnot(nrow(black) == 1L, nrow(white) == 1L)
    dp <- white$pseudobulk_log2cpm1 - black$pseudobulk_log2cpm1
    dm <- white$mean_lognorm - black$mean_lognorm
    reference <- pair_reference[pair_reference$population == definition$id & pair_reference$donor_group == donor_group, ]
    stopifnot(nrow(reference) == 1L, abs(dp - reference$diff_pseudobulk_log2cpm1) < 1e-12,
              abs(dm - reference$diff_mean_lognorm) < 1e-12,
              (sign(dp) != sign(dm)) == reference$metric_disagreement)
    data.frame(population = definition$id, donor_group, diff_pseudobulk_log2cpm1 = dp,
               diff_mean_lognorm = dm, metric_disagreement = sign(dp) != sign(dm))
  }))
}))
upper_limit <- ceiling(max(cells$rab7a_lognorm))
colors <- c(black = "#1a1a2e", white = "#c8c8c8")
plots <- lapply(prepared, function(item) {
  n_cells <- table(item$values$phenotype)
  labels <- c(black = paste0("black\n", format(n_cells[["black"]], big.mark = ","), " cells"),
              white = paste0("white\n", format(n_cells[["white"]], big.mark = ","), " cells"))
  plot <- ggplot(item$values, aes(phenotype, rab7a_lognorm, fill = phenotype)) +
    geom_violin(trim = TRUE, scale = "width", alpha = 0.7, linewidth = 0.4) +
    geom_point(position = position_jitter(width = 0.075, height = 0, seed = 20261001),
               size = 0.18, alpha = 0.12, color = "black") +
    geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.4, fill = "white", alpha = 0.95) +
    scale_fill_manual(values = colors) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(limits = c(0, upper_limit), expand = expansion(mult = c(0.015, 0.025))) +
    labs(title = item$definition$title, subtitle = "External cohort", x = NULL,
         y = "RAB7A log-normalized expression") +
    theme_classic(base_size = 10, base_family = font_family) +
    theme(legend.position = "none", axis.text = element_text(size = 9, color = "black"),
          axis.text.x = element_text(lineheight = 1.2),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 8.5, hjust = 0.5, margin = margin(b = 7)),
          plot.margin = margin(7, 8, 5, 5)) +
    annotate("text", x = 1.5, y = upper_limit * 0.96, label = item$p_label,
             family = font_family, size = 2.65, lineheight = 1.1, vjust = 0.5)
  built <- ggplot_build(plot)
  plotted_n <- vapply(split(built$data[[1]]$n, built$data[[1]]$group), function(n) unique(n)[1], numeric(1))
  stopifnot(identical(as.integer(plotted_n), as.integer(n_cells)),
            nrow(built$data[[2]]) == nrow(item$values), nrow(built$data[[3]]) == 2L,
            all(is.finite(built$data[[2]]$y)),
            identical(built$data[[2]]$y, item$values$rab7a_lognorm))
  plot
})
if (!dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output directory")
for (i in seq_along(prepared)) {
  base <- file.path(output_dir, prepared[[i]]$definition$filename)
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = 3.3, height = 4.1, family = font_family)
  print(plots[[i]])
  invisible(grDevices::dev.off())
  grDevices::png(paste0(base, ".png"), width = 3.3, height = 4.1, units = "in", res = 300,
                type = "cairo", bg = "white")
  print(plots[[i]])
  invisible(grDevices::dev.off())
}
write_tsv <- function(value, filename) {
  serialized <- lapply(value, function(column) {
    if (is.numeric(column)) sprintf("%.17g", column) else as.character(column)
  })
  write.table(as.data.frame(serialized), file.path(output_dir, filename), sep = "\t",
              quote = FALSE, row.names = FALSE)
}
statistics <- do.call(rbind, lapply(prepared, function(item) {
  data.frame(population = item$definition$id, n_black = item$definition$n[1], n_white = item$definition$n[2],
             statistic_w = unname(item$result$statistic), p_nominal = item$result$p.value,
             alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE,
             donor_adjusted = FALSE, multiple_testing_adjusted = FALSE)
}))
stopifnot(identical(input_hashes, tools::md5sum(input_paths)))
write_tsv(statistics, "recomputed_nominal_cell_tests.tsv")
write_tsv(paired, "recomputed_paired_descriptive.tsv")
write_tsv(data.frame(
  check = c("package_sha256_verified", "all_cells_and_zeros_plotted", "input_tables_unchanged", "nominal_statistics_match",
            "sample_summaries_match", "paired_summaries_match", "max_sample_absolute_difference",
            "jitter_seed", "font_family", "R_version", "ggplot2_version"),
  value = c("TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", sprintf("%.17g", max(sample_differences)),
            "20261001", font_family, as.character(getRversion()), as.character(packageVersion("ggplot2")))
), "replot_checks.tsv")
cat("Both panels reproduced with all cells and zeros. Nominal pooled-cell tests, not donor-adjusted.\n")
print(statistics[, c("population", "n_black", "n_white", "statistic_w", "p_nominal")], digits = 17)
