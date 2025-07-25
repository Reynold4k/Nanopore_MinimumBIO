# R Script for Processing Differential Coverage Data

# This script processes differential coverage data from experimental rounds and generates two types of plots: 
# a line plot depicting CPM values across experimental rounds
# a volcano plot highlighting significant changes in gene expression between conditions.

# What you what to modify:
# 1.Change the line 12 EXPERIMENTAL_FOLDER to your/exp/fastq parent path


# Set the experimental folder path
EXPERIMENTAL_FOLDER <- "D:/241004_exp/"

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)

# Split the path into components
path_components <- strsplit(EXPERIMENTAL_FOLDER, "/")[[1]]

# Extract the last three components
exp_name <- paste(rev(path_components)[1:3], collapse = "_")

# Get all round directories
round_dirs <- list.dirs(EXPERIMENTAL_FOLDER, full.names = TRUE, recursive = FALSE)

overall_CPM <- data.frame(Gene = character(), CPM = numeric(), Round = character(), stringsAsFactors = FALSE)

# Loop through each round directory
for (round_dir in round_dirs) {
  diff_coverage_file <- file.path(round_dir, "differential_coverage.txt")
  
  if (!file.exists(diff_coverage_file)) {
    warning(paste("Warning: File does not exist:", diff_coverage_file))
    next  # Skip to next iteration if the file does not exist
  }
  
  # Read the differential coverage file without header
  df <- read.table(diff_coverage_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Print the number of rows read and the first few lines to check the content
  print(paste("Number of rows read from", diff_coverage_file, ":", nrow(df)))
  print(head(df))  # Print the first few rows to check file content
  
  # Check if df has enough rows to set column names
  if (nrow(df) < 2) {
    warning(paste("Warning: Not enough data in", diff_coverage_file, ". Skipping."))
    next  # If there's not enough data, skip
  }
  
  # Set the last row as column names
  colnames(df) <- as.character(unlist(df[nrow(df), ]))
  
  # Remove the last row from the data as it contains the column names
  df <- df[-nrow(df), ]
  
  # Check if df is empty after processing
  if (nrow(df) == 0) {
    warning(paste("Warning: Empty dataframe for:", round_dir))
    next  # Skip this iteration if df is empty
  }
  
  # Convert columns to numeric as needed, but keep the first column (gene names) as character
  df <- df %>%
    mutate(across(-1, as.numeric))  # Convert all columns except the first one
  
  # Check the structure of the dataframe
  print(str(df))  # Print dataframe structure to check column names and data types
  
  # Calculate totals for control and experimental
  total_control <- sum(df$Control_Coverage, na.rm = TRUE)
  total_exp <- sum(df$Experimental_Coverage, na.rm = TRUE)
  
  # Calculate CPM_Control and CPM_EXP
  df <- df %>%
    mutate(CPM_Control = (Control_Coverage / total_control) * 1000000,
           CPM_EXP = (Experimental_Coverage / total_exp) * 1000000,
           Round = basename(round_dir))  # Add round information to the dataframe
  
  # Check if appropriate columns are present before proceeding
  if (!"Gene" %in% names(df) || !"CPM_EXP" %in% names(df)) {
    warning("Warning: Necessary columns 'Gene' or 'CPM_EXP' are missing in the processed dataframe.")
    next  # Skip to next round if necessary columns are missing
  }
  
  # Prepare data for overall CPM line plot
  overall_CPM <- bind_rows(overall_CPM,
                           data.frame(Gene = df$Gene, CPM = df$CPM_EXP - df$CPM_Control, Round = basename(round_dir))  # Assuming first column is 'Gene'
  )
  
  # Save updated dataframe
  updated_file_path <- file.path(round_dir, paste0("updated_", basename(round_dir), "_differential_coverage.txt"))
  write.table(df, updated_file_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  message(paste("Updated differential coverage saved to:", updated_file_path))
}

# Check for unique values in the Round column
unique_rounds <- unique(overall_CPM$Round)

# If there is only one unique round
if (length(unique_rounds) == 1) {
  original_round <- unique_rounds[1]  # Get the only round
  
  # Modify the round by subtracting one from the second character
  modified_round <- paste0(
    substr(original_round, 1, 1),  # "R"
    as.character(as.integer(substr(original_round, 2, 2)) - 1)  # Subtract 1 from the number
  )
  
  # Create new rows for each existing row with CPM set to 0
  new_rows <- overall_CPM %>%
    mutate(CPM = 0, Round = modified_round)
  
  # Append the new rows to overall_CPM
  overall_CPM <- bind_rows(overall_CPM, new_rows)
}

# Select top 11 most variable genes based on absolute mean CPM
top_genes <- overall_CPM %>%
  group_by(Gene) %>%
  summarise(mean_CPM = mean(CPM, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(abs(mean_CPM))) %>%
  slice_head(n = 11)  # Select top 11 based on the absolute value

# Filter overall_CPM for top genes to plot
top_genes_data <- overall_CPM %>%
  filter(Gene %in% top_genes$Gene)

# Reorder the Gene levels based on mean CPM for plotting
top_genes_data$Gene <- factor(top_genes_data$Gene, levels = top_genes$Gene[order(top_genes$mean_CPM, decreasing = TRUE)])

# Create the line plot only for the top 15 genes using CPM, with the legend ordered by mean CPM
line_plot <- ggplot(top_genes_data, aes(x = Round, y = CPM, group = Gene, color = Gene)) +
  geom_line() +
  geom_point() +
  xlab("Round") +
  ylab("CPM") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(
    plot.title = element_text(size = 14, face = "bold"), # Title size adjustment and bold
    legend.position = "right",  # Adjust as needed
    axis.title = element_text(size = 20),  # Axis title size
    axis.text = element_text(size = 18),  # Axis text size
    legend.title = element_text(size = 18),  # Legend title size
    legend.text = element_text(size = 14)   # Legend text size
  )

# Calculate log2foldchange and log10CPM for the last round's data
df <- df %>%
  mutate(
    log2foldchange = log2(CPM_EXP / CPM_Control),
    log10CPM = log10(abs(CPM_EXP - CPM_Control) + 1)  # 添加 +1 避免 log10(0)
  ) %>%
  # Delete the duplicated
  filter(!is.infinite(log2foldchange) & !is.infinite(log10CPM)) %>%
  mutate(
    color = case_when(
      log2foldchange < -1 ~ "blue",
      log2foldchange > 1 ~ "red",
      TRUE ~ "black"
    )
  )

# Calculate log2foldchange and log10CPM for the last round's data
df <- df %>%
  mutate(
    log2foldchange = log2(CPM_EXP / CPM_Control),
    log10CPM = log10(abs(CPM_EXP - CPM_Control) + 1)  # 添加 +1 避免 log10(0)
  ) %>%
  # Delete the duplicated
  filter(!is.infinite(log2foldchange) & !is.infinite(log10CPM)) %>%
  mutate(
    color = case_when(
      log2foldchange < -1 ~ "blue",
      log2foldchange > 1 ~ "red",
      TRUE ~ "black"
    )
  )

df <- df %>%
  filter(log10CPM >= 3)

# Select top 10 genes to label by log10CPM where the color is either red or blue
label_df <- df %>%
  filter(color %in% c("red", "blue")) %>%
  arrange(desc(log10CPM)) %>%
  slice_head(n = 100)  # Top 20 by log10CPM among significant changes


volcano_plot <- ggplot(df, aes(x = log2foldchange, y = log10CPM)) +
  geom_point(aes(color = color), alpha = 0.5) +
  scale_color_identity() +  # Directly use specified colors
  geom_text(data = label_df, aes(label = Gene),
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

ggsave(file.path(EXPERIMENTAL_FOLDER, "line_plot.png"), plot = line_plot, width = 8, height = 6)
ggsave(file.path(EXPERIMENTAL_FOLDER, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)
