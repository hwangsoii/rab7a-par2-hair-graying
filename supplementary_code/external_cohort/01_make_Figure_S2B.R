# Reproduce Figure S2B from the exact normalized source-data values released
# with this repository. No third-party FASTQ or full expression object is
# redistributed.

suppressPackageStartupMessages(library(ggplot2))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1) {
  normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
} else {
  normalizePath("01_make_Figure_S2B.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."),
                           mustWork = TRUE)
data_dir <- file.path(repo_root, "data", "external")
figure_dir <- Sys.getenv("HAIR_FIGURE_DIR",
                         unset = file.path(repo_root, "figures"))
result_dir <- Sys.getenv("HAIR_RESULTS_DIR",
                         unset = file.path(repo_root, "results"))
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

source_path <- file.path(data_dir, "Figure_S2B_source_data.tsv")
score_path <- file.path(data_dir, "Figure_S2B_cluster_scores.tsv")
sample_summary_path <- file.path(data_dir, "Figure_S2B_sample_summary.tsv")
stopifnot(file.exists(source_path), file.exists(score_path),
          file.exists(sample_summary_path))

source_data <- read.delim(source_path, stringsAsFactors = FALSE)
cluster_scores <- read.delim(score_path, stringsAsFactors = FALSE)
released_sample_summary <- read.delim(sample_summary_path,
                                      stringsAsFactors = FALSE)

required_source <- c("public_cell_id", "source_sample", "donor", "phenotype",
                     "krt_cluster", "KRT35", "KRT85",
                     "KRT35_KRT85_score", "RAB7A")
stopifnot(all(required_source %in% names(source_data)))
stopifnot(all(c("krt_cluster", "mean_krt35_85", "rank", "selected") %in%
                names(cluster_scores)))
stopifnot(!anyDuplicated(source_data$public_cell_id))
stopifnot(!anyNA(source_data[, required_source]))

expected_samples <- c("ryg029", "ryg035", "ryg047", "ryg048",
                      "black1_F62B", "white2_F62W")
stopifnot(setequal(unique(source_data$source_sample), expected_samples))
stopifnot(all(source_data$phenotype[source_data$source_sample == "ryg029"] ==
                "black"))

selected_clusters <- as.character(
  cluster_scores$krt_cluster[as.logical(cluster_scores$selected)]
)
expected_n_selected <- ceiling(nrow(cluster_scores) / 4)
stopifnot(length(selected_clusters) == expected_n_selected)
stopifnot(setequal(unique(as.character(source_data$krt_cluster)),
                   selected_clusters))

source_data$Phenotype <- factor(source_data$phenotype,
                                levels = c("black", "white"),
                                labels = c("Black", "White"))
black_values <- source_data$RAB7A[source_data$Phenotype == "Black"]
white_values <- source_data$RAB7A[source_data$Phenotype == "White"]
test_result <- wilcox.test(black_values, white_values, exact = FALSE)

statistics <- data.frame(
  comparison = "Black versus white",
  unit = "cell",
  n_black = length(black_values),
  n_white = length(white_values),
  mean_black = mean(black_values),
  mean_white = mean(white_values),
  median_black = median(black_values),
  median_white = median(white_values),
  wilcoxon_p = unname(test_result$p.value)
)
write.table(statistics,
            file.path(result_dir, "Figure_S2B_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sample_summary <- aggregate(
  RAB7A ~ source_sample + donor + phenotype,
  source_data,
  function(values) c(
    n_cells = length(values),
    mean = mean(values),
    median = median(values),
    pct_expressing = 100 * mean(values > 0)
  )
)
sample_summary <- data.frame(
  source_sample = sample_summary$source_sample,
  donor = sample_summary$donor,
  phenotype = sample_summary$phenotype,
  n_cells = as.integer(sample_summary$RAB7A[, "n_cells"]),
  mean_RAB7A = sample_summary$RAB7A[, "mean"],
  median_RAB7A = sample_summary$RAB7A[, "median"],
  pct_expressing_RAB7A = sample_summary$RAB7A[, "pct_expressing"]
)
stopifnot(isTRUE(all.equal(sample_summary, released_sample_summary,
                           check.attributes = FALSE, tolerance = 1e-10)))

p_label <- paste0("p = ", format.pval(test_result$p.value, digits = 2,
                                       eps = 1e-4))
plot_s2b <- ggplot(source_data, aes(x = Phenotype, y = RAB7A,
                                   fill = Phenotype)) +
  geom_violin(trim = TRUE, color = "black", linewidth = 0.35) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
               linewidth = 0.35) +
  geom_jitter(width = 0.13, size = 0.25, alpha = 0.15) +
  scale_fill_manual(values = c(Black = "#353544", White = "#E8E8E8"),
                    guide = "none") +
  coord_cartesian(ylim = c(0, max(source_data$RAB7A) * 1.13), clip = "on") +
  labs(title = "KRT35/KRT85-enriched keratinocytes", subtitle = p_label,
       x = NULL,
       y = "RAB7A expression") +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(figure_dir, "Figure_S2B.pdf"), plot_s2b,
       width = 3.6, height = 4.6, units = "in", device = cairo_pdf)

cat("Figure S2B written to:", file.path(figure_dir, "Figure_S2B.pdf"), "\n")
cat("Selected clusters:", paste(selected_clusters, collapse = ", "), "\n")
print(statistics)
