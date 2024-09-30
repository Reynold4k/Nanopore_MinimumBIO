#Installation of packages:
install.packages("dplyr")
install.packages("ggplot2")
install.packages("readr")


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)


# Set the experimental folder path
EXPERIMENTAL_FOLDER <- "/mnt/d/Bait_Glue/VHL/MB012/TON/230827"

# Split the path into components
path_components <- strsplit(EXPERIMENTAL_FOLDER, "/")[[1]]

# Extract the last three components
exp_name <- paste(rev(path_components)[1:3], collapse = "_")

# Get all round directories
round_dirs <- list.dirs(EXPERIMENTAL_FOLDER, full.names = TRUE, recursive = FALSE)

# Initialize variable to store overall CPM for line plotting
overall_CPM <- data.frame(Gene = character(), CPM = numeric(), Round = character(), stringsAsFactors = FALSE)

# Loop through each round directory
for (round_dir in round_dirs) {
  diff_coverage_file <- file.path(round_dir, "differential_coverage.txt")
  
  if (!file.exists(diff_coverage_file)) {
    warning(paste("Warning: File does not exist:", diff_coverage_file))
    next  # Skip to next iteration if the file does not exist
  }
  
  # Read the differential coverage file without header
  df <- read_table(diff_coverage_file, col_names = FALSE)
  
  # Set the last row as column names
  colnames(df) <- as.character(unlist(df[nrow(df), ]))
  
  # Remove the last row from the data as it contains the column names
  df <- df[-nrow(df), ]
  
  # Convert columns to numeric as needed, but keep the first column (gene names) as character
  df <- df %>%
    mutate(across(-1, as.numeric))  # Convert all columns except the first one
  
  # Calculate totals for control and experimental
  total_control <- sum(df$Control_Coverage, na.rm = TRUE)
  total_exp <- sum(df$Experimental_Coverage, na.rm = TRUE)
  
  # Calculate CPM_Control and CPM_EXP
  df <- df %>%
    mutate(CPM_Control = (Control_Coverage / total_control) * 1000000,
           CPM_EXP = (Experimental_Coverage / total_exp) * 1000000,
           Round = basename(round_dir))  # Add round information to the dataframe
  
  # Prepare data for overall CPM line plot
  overall_CPM <- bind_rows(overall_CPM,
                           data.frame(Gene = df$Gene, CPM = df$CPM_EXP, Round = basename(round_dir)))  # Use CPM_EXP for plotting
  
  # Save updated dataframe
  updated_file_path <- file.path(round_dir, paste0("updated_", basename(round_dir), "_differential_coverage.txt"))
  write.table(df, updated_file_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  message(paste("Updated differential coverage saved to:", updated_file_path))
}

# Select top 15 most variable genes based on absolute mean CPM
top_genes <- overall_CPM %>%
  group_by(Gene) %>%
  summarise(mean_CPM = mean(CPM, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(abs(mean_CPM))) %>%
  slice_head(n = 15)  # Select top 10 based on the absolute value

# Filter overall_CPM for top genes to plot
top_genes_data <- overall_CPM %>%
  filter(Gene %in% top_genes$Gene)

# Create the line plot only for the top 10 genes using CPM
line_plot <- ggplot(top_genes_data, aes(x = Round, y = CPM, group = Gene, color = Gene)) +
  geom_line() +
  geom_point() +
  labs(
    title = paste(exp_name),
    x = "Round",
    y = "CPM",
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(
    plot.title = element_text(size = 14, face = "bold"), # Title size adjustment and bold
    legend.position = "right",  # Adjust as needed
    axis.title = element_text(size = 20),  # Axis title size
    axis.text = element_text(size = 18),  # Axis text size
    legend.title = element_text(size = 18),  # Legend title size
    legend.text = element_text(size = 14)  # Legend text size
  )

ggsave(file.path(EXPERIMENTAL_FOLDER, "line_plot.png"), plot = line_plot, width = 8, height = 6)

# Calculate log2foldchange and log10CPM for the last round's data
df <- df %>%
  mutate(
    log2foldchange = log2(CPM_EXP / CPM_Control),
    log10CPM = log10(abs(CPM_EXP - CPM_Control) + 1), # Adding +1 to avoid log10(0)
    color = case_when(
      log2foldchange < -1 ~ "blue",
      log2foldchange > 1 ~ "red",
      TRUE ~ "black"
    )
  )

# Select top 10 genes to label by log10CPM where the color is either red or blue
label_df <- df %>%
  filter(color %in% c("red", "blue")) %>%
  arrange(desc(log10CPM)) %>%
  slice_head(n = 10)  # Top 20 by log10CPM among significant changes

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

ggsave(file.path(EXPERIMENTAL_FOLDER, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)
