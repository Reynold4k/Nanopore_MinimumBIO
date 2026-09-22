#!/usr/bin/env Rscript

# Consolidated differential gene-expression analysis for genome mode.
# Merges legacy Pipeline1/Pipeline1_limma.R and Pipeline1/Pipiline1_ttest.R
# with the documented v2 behaviour changes (C10, C16, C17, C18).

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
  library(readr)
})

# ---- CLI parser --------------------------------------------------------------

parse_args <- function(argv) {
  args <- list()
  i <- 1
  while (i <= length(argv)) {
    if (grepl("^--", argv[i])) {
      key <- sub("^--", "", argv[i])
      if (i + 1 <= length(argv) && !grepl("^--", argv[i + 1])) {
        args[[key]] <- argv[i + 1]
        i <- i + 2
      } else {
        args[[key]] <- TRUE
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  args
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

get_arg <- function(key, default = NULL, required = FALSE) {
  if (!key %in% names(args)) {
    if (required) stop("Missing required argument: --", key)
    return(default)
  }
  args[[key]]
}

exp_base_path      <- get_arg("exp", required = TRUE)
control_base_path  <- get_arg("control", required = TRUE)
gtf_file           <- get_arg("gtf", required = TRUE)
out_dir            <- get_arg("out", required = TRUE)
method             <- get_arg("method", default = "trend")
pvalue_threshold   <- as.numeric(get_arg("pvalue", default = "0.05"))
log10cpm_threshold <- as.numeric(get_arg("log10cpm", default = "3"))
top_label          <- as.integer(get_arg("top-label", default = "40"))
top_trajectory     <- as.integer(get_arg("top-trajectory", default = "10"))
round_pattern      <- get_arg("round-pattern", default = "R*")

if (!method %in% c("trend", "ttest")) {
  stop("--method must be either 'trend' or 'ttest'")
}

# ---- Input validation --------------------------------------------------------

if (!dir.exists(exp_base_path))      stop("Experiment directory does not exist: ", exp_base_path)
if (!dir.exists(control_base_path))  stop("Control directory does not exist: ", control_base_path)
if (!file.exists(gtf_file))          stop("GTF file does not exist: ", gtf_file)
if (!dir.exists(out_dir))            dir.create(out_dir, recursive = TRUE)

# ---- GTF annotation ----------------------------------------------------------

gtf_data <- import(gtf_file, format = "gtf")
genes <- gtf_data[gtf_data$type == "gene"]
gene_info <- data.frame(
  GeneID = as.character(mcols(genes)$gene_id),
  GeneName = as.character(mcols(genes)$gene_name),
  stringsAsFactors = FALSE
)
gene_info$GeneName <- ifelse(is.na(gene_info$GeneName) | gene_info$GeneName == "",
                             gene_info$GeneID, gene_info$GeneName)

# ---- Round discovery ---------------------------------------------------------

list_round_dirs <- function(base_path, pattern) {
  candidates <- Sys.glob(file.path(base_path, pattern))
  candidates[file.info(candidates)$isdir]
}

exp_dirs    <- list_round_dirs(exp_base_path, round_pattern)
control_dirs <- list_round_dirs(control_base_path, round_pattern)

common_dirs <- intersect(basename(exp_dirs), basename(control_dirs))
if (length(common_dirs) == 0) {
  stop("No common round directories found between --exp and --control matching pattern '",
       round_pattern, "'")
}

# ---- File paths --------------------------------------------------------------

generate_file_paths <- function(base_dir, common_dirs) {
  paths <- setNames(nm = common_dirs)
  for (dir in common_dirs) {
    file_path <- file.path(base_dir, dir, "step3", paste0(dir, "_combined_expression_counts.txt"))
    if (file.exists(file_path)) {
      paths[[dir]] <- file_path
    } else {
      warning("File does not exist: ", file_path)
    }
  }
  paths
}

control_paths <- generate_file_paths(control_base_path, common_dirs)
exp_paths     <- generate_file_paths(exp_base_path, common_dirs)

available <- common_dirs[!is.null(control_paths[common_dirs]) & !is.null(exp_paths[common_dirs])]
if (length(available) == 0) stop("No round with count files present in both experiment and control.")
control_paths <- control_paths[available]
exp_paths     <- exp_paths[available]
common_dirs   <- available

# ---- Read counts and compute CPM ---------------------------------------------

read_and_process <- function(paths) {
  cpm_results <- list()
  for (timepoint in names(paths)) {
    path <- paths[[timepoint]]
    counts <- read.delim(path, comment.char = "#", row.names = 1)
    n_annot <- min(5, ncol(counts))
    if (ncol(counts) <= 5) {
      stop("Count file has no sample columns after annotation columns: ", path)
    }
    counts_matrix <- as.matrix(counts[, -seq_len(n_annot), drop = FALSE])
    dge <- DGEList(counts = counts_matrix)
    cpm_results[[timepoint]] <- cpm(dge)
  }
  cpm_results
}

control_cpm_normalized <- read_and_process(control_paths)
exp_cpm_normalized     <- read_and_process(exp_paths)

# Align all CPM matrices to a common gene set across rounds
gene_sets <- c(lapply(control_cpm_normalized, rownames), lapply(exp_cpm_normalized, rownames))
common_genes <- Reduce(intersect, gene_sets)
if (length(common_genes) == 0) stop("No common genes found across rounds.")

control_cpm_normalized <- lapply(control_cpm_normalized, function(m) m[common_genes, , drop = FALSE])
exp_cpm_normalized     <- lapply(exp_cpm_normalized,     function(m) m[common_genes, , drop = FALSE])

# ---- Per-round exp - ctrl difference -----------------------------------------

calculate_difference <- function(experiment_cpm, control_cpm) {
  common <- intersect(rownames(experiment_cpm), rownames(control_cpm))
  exp_mean  <- rowMeans(experiment_cpm[common, , drop = FALSE])
  ctrl_mean <- rowMeans(control_cpm[common, , drop = FALSE])
  setNames(exp_mean - ctrl_mean, common)
}

differences <- sapply(seq_along(common_dirs), function(i) {
  calculate_difference(exp_cpm_normalized[[i]], control_cpm_normalized[[i]])
})
colnames(differences) <- common_dirs

differences_df <- as.data.frame(differences)
differences_df$Gene <- rownames(differences_df)

convert_to_gene_names <- function(differences_df, gene_info) {
  gene_names <- gene_info$GeneName[match(differences_df$Gene, gene_info$GeneID)]
  if (any(is.na(gene_names))) {
    warning("Some gene IDs could not be matched to gene names.")
  }
  differences_df$GeneName <- gene_names
  differences_df <- differences_df[!is.na(differences_df$GeneName), ]
  differences_df
}

differences_df_named <- convert_to_gene_names(differences_df, gene_info)

# ---- Variance and top genes for plotting -------------------------------------

variance_values <- apply(differences_df_named[, common_dirs, drop = FALSE], 1, sd, na.rm = TRUE)
differences_df_named$Variance <- variance_values

top_genes_line <- differences_df_named %>%
  arrange(desc(Variance)) %>%
  head(top_trajectory) %>%
  pull(GeneName)

long_df <- differences_df_named %>%
  pivot_longer(cols = all_of(common_dirs), names_to = "Timepoint", values_to = "Difference")

filtered_long_df <- long_df %>%
  filter(GeneName %in% top_genes_line)

filtered_long_df$Difference <- filtered_long_df$Difference / 1000

gene_order <- filtered_long_df %>%
  group_by(GeneName) %>%
  summarize(mean_expression = mean(Difference, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_expression))

filtered_long_df$GeneName <- factor(filtered_long_df$GeneName, levels = gene_order$GeneName)

colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
            "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
            "#8DA0CB", "#E5C494", "#B3B3B3", "#D84C3D", "#1F78B4",
            "#33A02C", "#FB9A99", "#A6CEE3", "#B2DF8A", "#FFEB3B")

path_components <- strsplit(exp_base_path, "/")[[1]]
exp_name <- paste(tail(path_components, 3), collapse = "_")

line_plot <- ggplot(filtered_long_df, aes(x = Timepoint, y = Difference, color = GeneName, group = GeneName)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = colors) +
  labs(
    title = exp_name,
    x = "Timepoints",
    y = "Normalized Gene Counts (Difference / 1000)",
    color = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)
  )

# ---- PCA (top 20 variable genes, legacy default) -----------------------------

topgenes_pca <- differences_df_named %>%
  arrange(desc(Variance)) %>%
  head(20) %>%
  pull(Gene)

get_avg_and_labels <- function(normalized_data, group_name, selected_genes) {
  timepoints <- names(normalized_data)
  avg_data_list <- list()
  labels_list <- list()
  for (timepoint in timepoints) {
    data <- normalized_data[[timepoint]]
    data <- data[rownames(data) %in% selected_genes, , drop = FALSE]
    avg_data <- rowMeans(data)
    avg_data_list[[timepoint]] <- avg_data
    labels_list[[timepoint]] <- data.frame(
      Sample = paste0(group_name, "_", timepoint),
      Timepoint = timepoint,
      Group = group_name,
      stringsAsFactors = FALSE
    )
  }
  list(data = do.call(rbind, avg_data_list), labels = do.call(rbind, labels_list))
}

control_avg <- get_avg_and_labels(control_cpm_normalized, "Control", topgenes_pca)
exp_avg     <- get_avg_and_labels(exp_cpm_normalized,     "Experiment", topgenes_pca)

combined_data  <- rbind(control_avg$data, exp_avg$data)
combined_labels <- rbind(control_avg$labels, exp_avg$labels)

pca_result <- prcomp(combined_data, center = TRUE, scale. = TRUE)
pca_data   <- cbind(as.data.frame(pca_result$x), combined_labels)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 5) +
  stat_ellipse(aes(group = Group, fill = Group),
               type = "norm",
               level = 0.95,
               alpha = 0.2,
               geom = "polygon") +
  scale_color_manual(values = c("Control" = "blue", "Experiment" = "red")) +
  scale_fill_manual(values = c("Control" = "blue", "Experiment" = "red")) +
  theme_classic() +
  labs(title = exp_name,
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14)
  )

# ---- Per-gene statistics (C16: mean over sample columns per round) -----------

time_points <- length(common_dirs)

exp_round_means  <- sapply(exp_cpm_normalized,     rowMeans)
ctrl_round_means <- sapply(control_cpm_normalized, rowMeans)
diff_round_means <- exp_round_means - ctrl_round_means

results <- data.frame(Gene = rownames(exp_cpm_normalized[[1]]), stringsAsFactors = FALSE)

if (method == "trend") {
  pvals  <- rep(NA_real_, nrow(results))
  slopes <- rep(NA_real_, nrow(results))
  t_vec  <- seq_len(time_points)
  for (i in seq_len(nrow(results))) {
    gene_id <- results$Gene[i]
    if (!(gene_id %in% rownames(diff_round_means))) next
    d <- diff_round_means[gene_id, ]
    fit <- lm(d ~ poly(t_vec, 1))
    coefs <- summary(fit)$coefficients
    if (nrow(coefs) >= 2) {
      pvals[i]  <- coefs[2, 4]
      slopes[i] <- coefs[2, 1]
    }
  }
  results$pvalue <- pvals
  results$slope  <- slopes
} else {
  pvals <- rep(NA_real_, nrow(results))
  for (i in seq_len(nrow(results))) {
    gene_id <- results$Gene[i]
    if (!(gene_id %in% rownames(exp_round_means))) next
    e <- exp_round_means[gene_id, ]
    c <- ctrl_round_means[gene_id, ]
    tt <- tryCatch(t.test(e, c, paired = FALSE), error = function(e) NULL)
    if (!is.null(tt)) pvals[i] <- tt$p.value
  }
  results$pvalue <- pvals
}

results$log_pvalue <- -log10(results$pvalue)

# Growth rate: mean exp CPM - mean ctrl CPM across all rounds
results$Growth_Rate <- rowMeans(exp_round_means) - rowMeans(ctrl_round_means)

results$normalized_Growth_Rate <- ifelse(
  results$Growth_Rate > 0,
  log10(results$Growth_Rate),
  -log10(-results$Growth_Rate + 1e-6)
)

# first/last-round control-subtracted average CPM (legacy growth-rate definition)
exp_last   <- exp_cpm_normalized[[length(exp_cpm_normalized)]]
exp_first  <- exp_cpm_normalized[[1]]
ctrl_last  <- control_cpm_normalized[[length(control_cpm_normalized)]]
ctrl_first <- control_cpm_normalized[[1]]

avg_last  <- rowMeans(exp_last)  - rowMeans(ctrl_last)
avg_first <- rowMeans(exp_first) - rowMeans(ctrl_first)

average_cpm <- avg_last - avg_first
results$log10_CPM    <- log10(average_cpm + 1e-6)
results$abs_log10_CPM <- abs(results$log10_CPM)

# BH-adjusted p-value (C10)
results$padj <- p.adjust(results$pvalue, method = "BH")

# Significance colouring at user-supplied p-value threshold (C10; not legacy limma 0.5 bug)
results$color_group <- ifelse(
  results$pvalue < pvalue_threshold,
  ifelse(results$Growth_Rate > 0, "Red", "Blue"),
  "Not Significant"
)

results$GeneName <- gene_info$GeneName[match(results$Gene, gene_info$GeneID)]

results <- results %>%
  filter(!is.na(normalized_Growth_Rate),
         !is.na(abs_log10_CPM),
         !is.na(color_group),
         !is.na(log_pvalue),
         abs_log10_CPM >= log10cpm_threshold)

# ---- Volcano plot ------------------------------------------------------------

highlight_genes <- results %>%
  filter(color_group %in% c("Red", "Blue"))

top_genes_CPM <- highlight_genes %>%
  arrange(desc(log10_CPM)) %>%
  slice_head(n = top_label)

top_genes_pvalue <- highlight_genes %>%
  arrange(desc(log_pvalue)) %>%
  slice_head(n = top_label)

top_genes_combined <- bind_rows(top_genes_CPM, top_genes_pvalue) %>%
  distinct()

top_genes_combined_filtered <- top_genes_combined[top_genes_combined$Gene %in% rownames(differences_df_named), ]

volcano_plot <- ggplot(results, aes(x = normalized_Growth_Rate, y = log_pvalue)) +
  geom_point(aes(color = color_group, size = abs_log10_CPM), alpha = 0.6) +
  scale_color_manual(
    values = c("Red" = "red", "Blue" = "blue", "Not Significant" = "grey"),
    breaks = c("Red", "Blue", "Not Significant"),
    labels = c("Significantly Upregulated", "Significantly Downregulated", "Not Significant")
  ) +
  theme_minimal() +
  labs(
    x = "Log10 Growth Rate", y = "-log10(p-value)",
    color = "Gene Regulation", size = "Absolute Log10 CPM"
  ) +
  geom_text_repel(
    data = top_genes_combined_filtered,
    aes(label = ifelse(!is.na(GeneName) & GeneName != "", GeneName, Gene)),
    size = 3, box.padding = 0.3, point.padding = 0.5,
    arrow = arrow(length = unit(0.01, "npc")),
    max.overlaps = Inf
  ) +
  scale_size_continuous(range = c(0, 6)) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 14),
    plot.margin = margin(5, 5, 5, 5, "pt")
  )

# ---- Save outputs ------------------------------------------------------------

differences_df_named$Gene <- rownames(differences_df_named)
rownames(differences_df_named) <- differences_df_named$Gene
differences_df_named <- differences_df_named %>%
  select(-Variance) %>%
  select(-Gene)

write.csv(differences_df_named, file.path(out_dir, "differences_df_named_no_variance.csv"), row.names = TRUE)

ggsave(file.path(out_dir, "pca_plot.png"),    plot = pca_plot,    width = 8,  height = 6)
ggsave(file.path(out_dir, "line_plot.png"),   plot = line_plot,   width = 8,  height = 6)
ggsave(file.path(out_dir, paste0(method, "_volcano_plot.png")), plot = volcano_plot, width = 12, height = 6)
