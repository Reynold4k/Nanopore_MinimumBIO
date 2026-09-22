#!/usr/bin/env Rscript

# Consolidated differential coverage analysis for targeted-reference mode.
# Merges legacy Pipeline1.5/Analysis.R and Analysis_withoutidmapping.R
# with documented v2 behaviour changes (C8, C17, C18, C20).

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

exp_base_path   <- get_arg("exp", required = TRUE)
out_dir         <- get_arg("out", required = TRUE)
id_mapping_file <- get_arg("id-mapping", default = NULL)

if (!dir.exists(exp_base_path)) stop("Experiment directory does not exist: ", exp_base_path)
if (!dir.exists(out_dir))       dir.create(out_dir, recursive = TRUE)

if (!is.null(id_mapping_file) && !file.exists(id_mapping_file)) {
  stop("ID mapping file does not exist: ", id_mapping_file)
}

# ---- Optional UniProt id mapping ---------------------------------------------

id_mapping <- NULL
if (!is.null(id_mapping_file)) {
  id_mapping <- read.table(id_mapping_file,
                           header = TRUE,
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           fill = TRUE,
                           quote = "",
                           comment.char = "")
  if (!all(c("From", "Gene.Names") %in% names(id_mapping))) {
    warning("ID mapping file does not contain required columns 'From' and 'Gene.Names'; mapping skipped.")
    id_mapping <- NULL
  }
}

apply_id_mapping <- function(df, id_map) {
  if (is.null(id_map) || nrow(df) == 0) return(df)
  df %>%
    mutate(Gene_Match = sub("^([^_]*_[^_]*)_.*$", "\\1", Gene)) %>%
    left_join(id_map, by = c("Gene_Match" = "From")) %>%
    mutate(Gene = ifelse(is.na(Gene.Names), Gene,
                         paste(Gene.Names, sub("^[^_]*_[^_]*_", "", Gene), sep = "_"))) %>%
    select(-Gene_Match, -Gene.Names)
}

# ---- Round loop --------------------------------------------------------------

path_components <- strsplit(exp_base_path, "/")[[1]]
exp_name <- paste(tail(path_components, 3), collapse = "_")

round_dirs <- list.dirs(exp_base_path, full.names = TRUE, recursive = FALSE)
if (length(round_dirs) == 0) stop("No round directories found under: ", exp_base_path)

overall_CPM <- data.frame(Gene = character(), CPM = numeric(), Round = character(), stringsAsFactors = FALSE)
last_df <- NULL

for (round_dir in round_dirs) {
  diff_coverage_file <- file.path(round_dir, "differential_coverage.txt")

  if (!file.exists(diff_coverage_file)) {
    warning("File does not exist: ", diff_coverage_file)
    next
  }

  # C8: v2 Compare writes a proper header; read with header=TRUE.
  df <- read.table(diff_coverage_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  message("Read ", nrow(df), " rows from ", diff_coverage_file)

  # Tolerate legacy trailing-header format: drop a last row equal to the column names.
  if (nrow(df) >= 1) {
    last_row <- as.character(unlist(df[nrow(df), ]))
    if (length(last_row) == ncol(df) && all(last_row == colnames(df))) {
      df <- df[-nrow(df), ]
      message("Removed trailing header row from ", basename(round_dir))
    }
  }

  if (nrow(df) == 0) {
    warning("Empty dataframe after header handling for: ", round_dir)
    next
  }

  if (!all(c("Gene", "Control_Coverage", "Experimental_Coverage") %in% names(df))) {
    warning("Required columns missing in ", diff_coverage_file, "; skipping.")
    next
  }

  df <- df %>%
    mutate(across(c(Control_Coverage, Experimental_Coverage), as.numeric))

  total_control <- sum(df$Control_Coverage, na.rm = TRUE)
  total_exp     <- sum(df$Experimental_Coverage, na.rm = TRUE)

  if (total_control == 0 || total_exp == 0) {
    warning("Zero total coverage in ", basename(round_dir), "; skipping.")
    next
  }

  df <- df %>%
    mutate(CPM_Control = (Control_Coverage / total_control) * 1000000,
           CPM_EXP     = (Experimental_Coverage / total_exp) * 1000000,
           Round       = basename(round_dir))

  df <- apply_id_mapping(df, id_mapping)

  overall_CPM <- bind_rows(
    overall_CPM,
    data.frame(Gene = df$Gene,
               CPM  = df$CPM_EXP - df$CPM_Control,
               Round = basename(round_dir),
               stringsAsFactors = FALSE)
  )

  updated_file_path <- file.path(out_dir, paste0("updated_", basename(round_dir), "_differential_coverage.txt"))
  write.table(df, updated_file_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  message("Updated differential coverage saved to: ", updated_file_path)

  last_df <- df
}

if (nrow(overall_CPM) == 0) stop("No valid differential coverage data found.")

overall_CPM <- apply_id_mapping(overall_CPM, id_mapping)

# ---- Line plot (C17: skip if fewer than 2 rounds) ----------------------------

n_rounds <- length(unique(overall_CPM$Round))

if (n_rounds < 2) {
  message("Skipping line plot: fewer than 2 rounds available (", n_rounds, " round).")
  line_plot <- NULL
} else {
  top_genes <- overall_CPM %>%
    group_by(Gene) %>%
    summarise(mean_CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(abs(mean_CPM))) %>%
    slice_head(n = 10)

  top_genes_data <- overall_CPM %>%
    filter(Gene %in% top_genes$Gene)

  top_genes_data$Gene <- factor(top_genes_data$Gene,
                                levels = top_genes$Gene[order(top_genes$mean_CPM, decreasing = TRUE)])

  line_plot <- ggplot(top_genes_data, aes(x = Round, y = CPM, group = Gene, color = Gene)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    xlab("Round") +
    ylab("CPM") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right",
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 14)
    )
}

# ---- Volcano plot (last round; C20: single mutate block) ---------------------

if (is.null(last_df) || nrow(last_df) == 0) {
  stop("No round data available for volcano plot.")
}

last_df <- last_df %>%
  mutate(
    log2foldchange = log2(CPM_EXP / CPM_Control),
    log10CPM       = log10(abs(CPM_EXP - CPM_Control) + 1)
  ) %>%
  filter(!is.infinite(log2foldchange) & !is.infinite(log10CPM)) %>%
  mutate(
    color = case_when(
      log2foldchange < -1 ~ "blue",
      log2foldchange >  1 ~ "red",
      TRUE                ~ "black"
    )
  )

top_log10CPM <- last_df %>%
  arrange(desc(abs(log10CPM))) %>%
  slice_head(n = 30)

top_log2foldchange <- last_df %>%
  arrange(desc(abs(log2foldchange))) %>%
  slice_head(n = 30)

label_df <- bind_rows(top_log10CPM, top_log2foldchange) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  filter(color != "black")

volcano_plot <- ggplot(last_df, aes(x = log2foldchange, y = log10CPM)) +
  geom_point(aes(color = color), alpha = 0.5, size = 4) +
  scale_color_identity() +
  geom_text(data = label_df,
            aes(label = Gene),
            size = 4, vjust = -0.5, hjust = 0.5, check_overlap = TRUE) +
  labs(title = paste("Volcano Plot of", exp_name),
       x = "Log2 Fold Change",
       y = "Log10 CPM") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# ---- Save outputs ------------------------------------------------------------

if (!is.null(line_plot)) {
  ggsave(file.path(out_dir, "line_plot.png"), plot = line_plot, width = 8, height = 6)
}
ggsave(file.path(out_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)
