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



id_mapping <- read.table("D:/Bait_Glue/2024_10_15.tsv", 
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
  print(head(df)) 
  
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
                           data.frame(Gene = df$Gene, CPM = df$CPM_EXP - df$CPM_Control, Round = basename(round_dir))  # Assuming first column is 'Gene'
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

overall_CPM <- overall_CPM %>%
  mutate(Gene_Match = sub("^([^_]*_[^_]*)_.*$", "\\1", Gene))  


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
  # 1. ggplot(...): Initializes a ggplot object for plotting.
  #    - top_genes_data: This specifies the data frame used for the plot.
  #    - aes(...): Aesthetic mappings: 'Round' is mapped to the x-axis and 'CPM' (Counts Per Million) is mapped to the y-axis.
  #      Each 'Gene' is represented as a different color and grouped by 'Gene'.

  geom_line(linewidth = 1.2) +  # Set line thickness  
  geom_point(size = 3) +  # Set point size

  geom_line() +
  # 2. geom_line(): Adds lines connecting data points for each 'Gene', tracking progression across 'Round'.
  
  geom_point() +
  # 3. geom_point(): Adds points at each data position on top of the lines, making individual data points visually distinct.
  
  xlab("Round") +
  # 4. xlab("Round"): Sets the label for the x-axis to "Round".
  
  ylab("CPM") +
  # 5. ylab("CPM"): Sets the label for the y-axis to "CPM", which stands for counts per million, a common normalization form in RNA-seq data.
  
  theme_minimal() +
  # 6. theme_minimal(): Applies a minimal theme to the plot, providing a clean, uncluttered aesthetic with simple grid lines.
  
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    # 7. axis.text.x = element_text(angle = 90, hjust = 1):
    #    - Rotates x-axis text by 90 degrees to improve readability when there are many levels.
    #    - hjust = 1 right-aligns the text, useful when text is rotated.
    
    plot.title = element_text(size = 14, face = "bold"),
    # 8. plot.title: Sets the size of the plot title text to 14 and makes it bold for emphasis.
    
    legend.position = "right",
    # 9. legend.position = "right": Positions the legend on the right side of the plot for easy access and readability.
    
    axis.title = element_text(size = 20),
    # 10. axis.title: Sets the text size of both x and y-axis titles to 20, enhancing their prominence.
    
    axis.text = element_text(size = 18),
    # 11. axis.text: Sets the text size for axis labels (tick marks) to 18 for ease of reading.
    
    legend.title = element_text(size = 18),
    # 12. legend.title: Configures the legend title text size to 18, ensuring it is easily readable.
    
    legend.text = element_text(size = 14),
    # 13. legend.text: Sets the size of the text in the legend entries to 14, maintaining clarity.
    
    # Set transparent backgrounds
    panel.background = element_rect(fill = NA, color = NA),
    # 14. panel.background: Removes the background fill and border from the plot's main panel, making it transparent.
    
    plot.background = element_rect(fill = NA, color = NA),
    # 15. plot.background: Removes the background fill and border from the entire plot area, enabling full transparency.
    
    legend.background = element_rect(fill = NA, color = NA)
    # 16. legend.background: Makes the background behind the legend transparent, removing any fill or border.
  )



df <- df %>%
  # Start with the existing data frame `df` and apply a series of transformations using the pipe operator `%>%`.
  
  mutate(
    log2foldchange = log2(CPM_EXP / CPM_Control),
    # 1. Calculate the log2 fold change:
    #    - CPM_EXP / CPM_Control: Computes the ratio of counts per million (CPM) between Experimental and Control conditions.
    #    - log2(...): Takes the base-2 logarithm of the ratio. This transformation is symmetrical for values >1 and <1, making it intuitive for interpreting up and down regulation.

    log10CPM = log10(abs(CPM_EXP - CPM_Control) + 1)
    # 2. Calculate log10 CPM difference:
    #    - abs(CPM_EXP - CPM_Control): Computes the absolute difference between CPM values of Experimental and Control to avoid negative values.
    #    - +1: Prevents taking the logarithm of zero, which is undefined. Add 1 to shift the scale.
    #    - log10(...): Takes the base-10 logarithm of the adjusted CPM difference, scaling the data for better visualization and interpretation.
  ) %>%
  
  filter(!is.infinite(log2foldchange) & !is.infinite(log10CPM)) %>%
  # 3. Remove infinite values:
  #    - filter(...): Excludes any rows where `log2foldchange` or `log10CPM` are infinite.
  #    - This typically occurs when dividing by zero or encountering missing values, ensuring that subsequent analyses are performed only on finite, meaningful values.

  mutate(
    color = case_when(
      log2foldchange < -1 ~ "blue",
      # 4. Assign colors based on fold change:
      #    - log2foldchange < -1: Assigns "blue" to indicate significant down-regulation, where the expression is higher in the Control than in the Experimental condition.

      log2foldchange > 1 ~ "red",
      #    - log2foldchange > 1: Assigns "red" to represent significant up-regulation, where the expression is higher in the Experimental than in the Control condition.

      TRUE ~ "black"
      #    - TRUE: Assigns "black" when conditions above aren't met, indicating no significant regulation (change is between -1 and 1 in log2 fold change).
    )
  )

# Select top 30 genes by absolute log10CPM
top_log10CPM <- df %>%
  arrange(desc(abs(log10CPM))) %>%
  slice_head(n = 30)

# Select top 30 genes by absolute log2foldchange
top_log2foldchange <- df %>%
  arrange(desc(abs(log2foldchange))) %>%
  slice_head(n = 30)

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
