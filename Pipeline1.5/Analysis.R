# R Script for Processing Differential Coverage Data

# This script processes differential coverage data from experimental rounds and generates two types of plots: 
# a line plot depicting CPM values across experimental rounds
# a volcano plot highlighting significant changes in gene expression between conditions.

# What you what to modify:
# 1.Change the line 19 EXPERIMENTAL_FOLDER to your/exp/fastq parent path
# 2.Change the line 23 id_mapping to your/id_mapping_file


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)


# Set the experimental folder path
EXPERIMENTAL_FOLDER <- "/srv/scratch/z3546698/tutorial/Small_Molecule/FK506/T7pep/240824"



id_mapping <- read.table("/srv/scratch/z3546698/tutorial/reference/idmapping_2024_10_01.tsv", 
                         header = TRUE, 
                         sep = "\t", 
                         stringsAsFactors = FALSE, 
                         fill = TRUE, 
                         quote = "",  
                         comment.char = "") 





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
  print(head(df))  # 打印前几行查看文件内容
  
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
  print(str(df))  
  
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
                           data.frame(Gene = df$Gene, CPM = df$CPM_EXP, Round = basename(round_dir))  # Assuming first column is 'Gene'
  )
  
  # Save updated dataframe
  updated_file_path <- file.path(round_dir, paste0("updated_", basename(round_dir), "_differential_coverage.txt"))
  write.table(df, updated_file_path, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  message(paste("Updated differential coverage saved to:", updated_file_path))
}



df <- df %>%
  mutate(Gene_Match = sub("^([^_]*_[^_]*)_.*$", "\\1\\2", Gene))  

# left connection with id_mapping
df <- df %>%
  left_join(id_mapping, by = c("Gene_Match" = "From")) %>%
  mutate(Gene = ifelse(is.na(Gene.Names), Gene, 
                       paste(Gene.Names, sub("^[^_]*_[^_]*_", "", Gene), sep="_"))) 

# 提取第二个下划线前的内容
overall_CPM <- overall_CPM %>%
  mutate(Gene_Match = sub("^([^_]*_[^_]*)_.*$", "\\1", Gene))  # 提取到第二个下划线前的内容

# 使用 Gene_Match 列与 id_mapping 进行左连接
overall_CPM <- overall_CPM %>%
  left_join(id_mapping, by = c("Gene_Match" = "From")) %>%
  mutate(Gene = ifelse(is.na(Gene.Names), Gene, 
                       paste(Gene.Names, sub("^[^_]*_[^_]*_", "", Gene), sep="_")))

# Select top 15 most variable genes based on absolute mean CPM
top_genes <- overall_CPM %>%
  group_by(Gene) %>%
  summarise(mean_CPM = mean(CPM, na.rm = TRUE), .groups = 'drop') %>%
  arrange(desc(abs(mean_CPM))) %>%
  slice_head(n = 10)  # Select top 15 based on the absolute value

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

# Select top 10 genes by absolute log10CPM
top_log10CPM <- df %>%
  arrange(desc(abs(log10CPM))) %>%
  slice_head(n = 10)

# Select top 10 genes by absolute log2foldchange
top_log2foldchange <- df %>%
  arrange(desc(abs(log2foldchange))) %>%
  slice_head(n = 10)

# Combine both lists and remove duplicates
label_df <- bind_rows(top_log10CPM, top_log2foldchange) %>%
  distinct(Gene, .keep_all = TRUE)

# Create the volcano plot
volcano_plot <- ggplot(df, aes(x = log2foldchange, y = log10CPM)) +
  # 1. ggplot(): Initializes a ggplot object for plotting.
  #    - df: The dataframe containing the data to be plotted.
  #    - aes(...): Aesthetic mappings, defining 'log2foldchange' as the x variable and 'log10CPM' as the y variable for the plot's axes.
  
  geom_point(aes(color = color), alpha = 0.5, size = 5) +
  # 2. geom_point(...): Adds a layer for plotting points.
  #    - color = color: Maps the 'color' column in the data to the color of the points.
  #    - alpha = 0.5: Sets the transparency level of the points (0 is fully transparent, 1 is fully opaque).
  #    - size = 5: Increases the size of each point, using ggplot's size scale.
  
  scale_color_identity() +
  # 3. scale_color_identity(): Uses the exact colors specified in the 'color' aesthetic without scaling.
  
  geom_text(data = label_df, aes(label = Gene),
            size = 4, vjust = -0.5, hjust = 0.5, check_overlap = TRUE) +
  # 4. geom_text(...): Adds text labels to each point.
  #    - data = label_df: Uses a separate dataframe for the labels, typically a subset of points.
  #    - aes(label = Gene): Specifies that the text for each label should come from the 'Gene' column.
  #    - size = 4: Sets the size of the text labels.
  #    - vjust = -0.5: Vertically adjusts text to be slightly above the point.
  #    - hjust = 0.5: Horizontally centers the text on the point.
  #    - check_overlap = TRUE: Ensures that text labels do not overlap, omits some labels if needed.
  
  labs(title = paste("Volcano Plot of", exp_name),
       x = "Log_2 Fold Change",
       y = "Log_10 CPM") +
  # 5. labs(...): Modifies the plot labels such as the title and axis labels.
  #    - title: The plot title, dynamically combines "Volcano Plot of" with the variable 'exp_name'.
  #    - x: The label for the x-axis.
  #    - y: The label for the y-axis.
  
  ylim(1, 6) +
  # 6. ylim(1, 6): Sets the limits for the y-axis between 1 and 6.
  #    - Values outside this range will not be plotted on the y-axis.
  
  theme_minimal() +
  # 7. theme_minimal(): Applies a minimalistic theme with clean, simple grid lines and no background shading.
  
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
# 8. theme(...): Customizes text and theme elements.
#    - plot.title: Sets the style of the plot title, with size 14 and bold text.
#    - axis.title: Sets the size of the axis titles, size 18.
#    - axis.text: Sets the size of the axis tick labels, making them very large at size 24.
#    - legend.title: Sets the size of the legend title text, size 14.
#    - legend.text: Sets the size of the legend item text, size 12.




ggsave(file.path(EXPERIMENTAL_FOLDER, "line_plot.png"), plot = line_plot, width = 8, height = 6)
ggsave(file.path(EXPERIMENTAL_FOLDER, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)
