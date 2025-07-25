
# Define paths
exp_base_path <- "D:/A_FKBP1B/MB022/CoT/241018"
control_base_path <- "D:/A_FKBP1B/No_glue/CoT/241018"


# Read GTF file
gtf_file <- "D:/hg38/Homo_sapiens.GRCh38.112.gtf"  # Please replace with the actual path


library(ggrepel)
library(Matrix)
library(DelayedArray)
library(MASS)
library(mgcv)
library(SummarizedExperiment)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(dplyr)
library(tidyr)

#Path to your newly generated Routput folder


# Output paths
plot_base_dir <- file.path(exp_base_path, "plots")  # Construct path to Routput in the experiment directory


gtf_data <- import(gtf_file, format="gtf")

# Extract gene information
genes <- gtf_data[gtf_data$type == "gene"]

# Extract gene ID and gene name
gene_info <- data.frame(
  GeneID = sapply(mcols(genes)$gene_id, function(x) x),
  GeneName = sapply(mcols(genes)$gene_name, function(x) x),
  stringsAsFactors = FALSE
)

# Convert gene IDs to gene names
convert_to_gene_names <- function(cpm_matrix, gene_info) {
  gene_ids <- rownames(cpm_matrix)
  cpm_matrix$GeneName <- gene_info$GeneName[match(gene_ids, gene_info$GeneID)]
  rownames(cpm_matrix) <- cpm_matrix$GeneName
  cpm_matrix <- cpm_matrix[, -which(names(cpm_matrix) == "GeneName")]
  return(cpm_matrix)
}


# Load required libraries
library(edgeR)

# Define function to list directories
list_dirs <- function(base_path) {
  dirs <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[basename(dirs) %in% paste0("R", 0:9)]  
  return(dirs)
}


# Split the path into components
path_components <- strsplit(exp_base_path, "/")[[1]]

# Extract the last three components
exp_name <- paste(rev(path_components)[1:3], collapse = "_")

# Get list of directories for each path
exp_dirs <- list_dirs(exp_base_path)
control_dirs <- list_dirs(control_base_path)

# Determine the number of replicates (timepoints) to use
common_dirs <- intersect(basename(exp_dirs), basename(control_dirs))
max_timepoints <- length(common_dirs)

# Generate the file paths for each replicate that exists in both experiment and control
generate_file_paths <- function(base_dir, common_dirs) {
  paths <- setNames(nm = common_dirs)
  for (dir in common_dirs) {
    # Construct the path to the step4 directory within each timepoint
    step3_dir <- file.path(base_dir, dir, "step3")
    
    # Construct the full file path to the expression counts file
    file_path <- file.path(step3_dir, paste0(dir, "_combined_expression_counts.txt"))
    
    # Check if the file exists and store the path, otherwise, issue a warning
    if (file.exists(file_path)) {
      paths[[dir]] <- file_path
    } else {
      warning(paste("File does not exist:", file_path))
    }
  }
  return(paths)
}

# Generate full paths
control_paths <- generate_file_paths(control_base_path, common_dirs)
exp_paths <- generate_file_paths(exp_base_path, common_dirs)

read_and_process <- function(paths) {
  cpm_results <- list()  # Initialize a list to store CPM normalized matrices for each timepoint
  
  for (timepoint in names(paths)) {
    path <- paths[[timepoint]]
    
    if (file.exists(path)) {
      # Read and process the count data
      counts <- read.delim(path, comment.char="#", row.names=1)
      counts_matrix <- as.matrix(counts[, -c(1:5)])  # Adjust column indices based on file format
      
      # Create DGEList object and calculate CPM
      dge <- DGEList(counts=counts_matrix)
      normalized_matrix <- cpm(dge)
      
      # Store the normalized matrix in the list with the key being the timepoint
      cpm_results[[timepoint]] <- normalized_matrix
    } else {
      warning(paste("File does not exist:", path))
    }
  }
  
  return(cpm_results)
}

# Assuming control_paths and exp_paths are defined as in previous sections
control_cpm_normalized <- read_and_process(control_paths)
exp_cpm_normalized <- read_and_process(exp_paths)


# Align and adjust data
calculate_difference <- function(experiment_cpm, control_cpm) {
  # Ensure gene names match in both matrices
  common_genes <- intersect(rownames(experiment_cpm), rownames(control_cpm))
  
  # Align matrices
  exp_matrix <- experiment_cpm[common_genes, , drop=FALSE]
  ctrl_matrix <- control_cpm[common_genes, , drop=FALSE]
  
  # Calculate the mean for each gene
  exp_mean <- rowMeans(exp_matrix)
  ctrl_mean <- rowMeans(ctrl_matrix)
  
  # Calculate difference
  difference <- exp_mean - ctrl_mean
  
  return(difference)
}

# Calculate difference for each time point
timepoints <- common_dirs

differences <- sapply(seq_along(timepoints), function(i) {
  calculate_difference(exp_cpm_normalized[[i]],
                       control_cpm_normalized[[i]])
})

colnames(differences) <- timepoints

# Convert matrix to data frame
differences_df <- as.data.frame(differences)
differences_df$Gene <- rownames(differences)



# Convert gene IDs to gene names
convert_to_gene_names <- function(differences_df, gene_info) {
  # Ensure both data frames have GeneID column
  gene_ids <- differences_df$Gene
  
  # Update gene information
  gene_names <- gene_info$GeneName[match(gene_ids, gene_info$GeneID)]
  
  # Ensure all gene IDs have corresponding names
  if (any(is.na(gene_names))) {
    warning("Some gene IDs could not be matched to gene names.")
  }
  
  # Add gene names to data frame
  differences_df$GeneName <- gene_names
  
  # Remove rows that did not match a name
  differences_df <- differences_df[!is.na(differences_df$GeneName), ]
  
  # Optionally remove the old Gene column
  # differences_df$Gene <- NULL
  
  return(differences_df)
}

# Use the updated function to convert
differences_df_named <- convert_to_gene_names(differences_df, gene_info)

# Convert data frame to long format
long_df <- differences_df_named %>%
  pivot_longer(cols = starts_with("R"), names_to = "Timepoint", values_to = "Difference")

# Calculate variance for each gene (e.g., standard deviation across timepoints)
calculate_variance <- function(differences_df) {
  variance_values <- apply(differences_df, 1, sd, na.rm = TRUE)
  return(variance_values)
}

# Calculate variance for each gene
variance_values <- calculate_variance(as.matrix(differences_df_named[, -ncol(differences_df_named)]))

# Add variance to data frame
differences_df_named$Variance <- variance_values

# Select the top 10 genes with the largest variance
top_genes <- differences_df_named %>%
  arrange(desc(Variance)) %>%
  head(10) %>%
  pull(GeneName)

# Filter data for top 20 genes
filtered_long_df <- long_df %>%
  filter(GeneName %in% top_genes)

filtered_long_df$Difference <- filtered_long_df$Difference / 1000

topgenes <- differences_df_named %>%
  arrange(desc(Variance)) %>%
  head(20) %>%
  pull(Gene)

# Calculate averages and create labels
get_avg_and_labels <- function(normalized_data, group_name) {
  timepoints <- names(normalized_data)
  avg_data_list <- list()
  labels_list <- list()
  
  for (timepoint in timepoints) {
    data <- normalized_data[[timepoint]]
    # Select genes of interest
    data <- data[rownames(data) %in% topgenes, ]
    
    # Compute average for each gene
    avg_data <- rowMeans(data)
    
    # Save averages and labels
    avg_data_list[[timepoint]] <- avg_data
    labels_list[[timepoint]] <- data.frame(
      Sample = paste0(group_name, "_", timepoint),
      Timepoint = timepoint,
      Group = group_name
    )
  }
  
  avg_data_matrix <- do.call(rbind, avg_data_list)
  labels_df <- do.call(rbind, labels_list)
  
  return(list(data = avg_data_matrix, labels = labels_df))
}

# Get average data and labels from control and exp
control_avg <- get_avg_and_labels(control_cpm_normalized, "Control")
exp_avg <- get_avg_and_labels(exp_cpm_normalized, "Experiment")

# Combine all data and labels
combined_data <- rbind(control_avg$data, exp_avg$data)
combined_labels <- rbind(control_avg$labels, exp_avg$labels)

pca_result <- prcomp(combined_data, center = TRUE, scale. = TRUE)

# Extract the PCA scores
pca_scores <- as.data.frame(pca_result$x)

# Combine the PCA scores with the labels
# combined_labels should contain the metadata necessary for the plot
pca_data <- cbind(pca_scores, combined_labels)



pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 5) +
  stat_ellipse(aes(group = Group, fill = Group), 
               type = "norm", 
               level = 0.95, 
               alpha = 0.2,  # Set transparency
               geom = "polygon") +  # Fill ellipse with polygon
  scale_color_manual(values = c("Control" = "blue", "Experiment" = "red")) +
  scale_fill_manual(values = c("Control" = "blue", "Experiment" = "red")) +  # Set fill color
  theme_classic() +
  labs(title = paste(exp_name),
       x = "Principal Component 1",
       y = "Principal Component 2")+
  theme(
    plot.title = element_text(size = 14, face = "bold"), # Title size adjustment and bold
    legend.position = "right",  # Adjust as needed
    axis.title = element_text(size = 20),  # Axis title size
    axis.text = element_text(size = 18),  # Axis text size
    legend.title = element_text(size = 18),  # Legend title size
    legend.text = element_text(size = 14)  # Legend text size
  )


library(ggplot2)

# Define 20 colors
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#FFFF33",
            "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
            "#8DA0CB", "#E5C494", "#B3B3B3", "#D84C3D", "#1F78B4",
            "#33A02C", "#FB9A99", "#A6CEE3", "#B2DF8A", "#FFEB3B")

# Draw line plot and set colors with scale_color_manual

# Compute the ordering metric (e.g., mean expression level) for each gene
gene_order <- filtered_long_df %>%
  group_by(GeneName) %>%
  summarize(mean_expression = mean(Difference, na.rm = TRUE)) %>%
  arrange(desc(mean_expression))

# Reorder the GeneName factor levels based on the mean expression
filtered_long_df$GeneName <- factor(filtered_long_df$GeneName, levels = gene_order$GeneName)

# Create the plot with the reordered legend
line_plot <- ggplot(filtered_long_df, aes(x = Timepoint, y = Difference, color = GeneName, group = GeneName)) +
  geom_line(size = 1.2) +  # Set line thickness
  geom_point(size = 3) +  # Set point size
  scale_color_manual(values = colors) +  # Use custom colors
  labs(
    title = paste(exp_name),
    x = "Timepoints",
    y = "Normalized Gene Counts",
    color = "Gene"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"), # Title size adjustment and bold
    legend.position = "right",  # Adjust as needed
    axis.title = element_text(size = 20),  # Axis title size
    axis.text = element_text(size = 18),  # Axis text size
    legend.title = element_text(size = 18),  # Legend title size
    legend.text = element_text(size = 14)  # Legend text size
  )

timepoints <- names(exp_cpm_normalized)  # Get a list of all time points from exp_cpm_normalized

time_points <- length(timepoints)  

# Create a result data frame to store p-values
results <- data.frame(Gene = rownames(exp_cpm_normalized[[1]]), pvalue = NA)

# Loop through each gene
for (gene_index in 1:nrow(exp_cpm_normalized[[1]])) {
  # Extract the expression values for this gene at all time points
  exp_values <- sapply(exp_cpm_normalized, function(x) x[gene_index, ])  # Extract from experimental data
  control_values <- sapply(control_cpm_normalized, function(x) x[gene_index, ])  # Extract from control data
  # Perform unpaired t-test
  t_test_result <- t.test(exp_values, control_values, paired = FALSE)
  # Extract and save p-value
  results$pvalue[gene_index] <- t_test_result$p.value
}



# Combine the results with the original gene list
results$Gene <- rownames(exp_cpm_normalized[[1]])

# Next, you can perform significance analysis, calculate growth rates, or generate volcano plots
# For example:
results$log_pvalue <- -log10(results$pvalue)

# Calculate growth rate
results$Growth_Rate <- rowMeans(simplify2array(exp_cpm_normalized)) - rowMeans(simplify2array(control_cpm_normalized))

# Log10 transformation for growth rate
results$normalized_Growth_Rate <- ifelse(
  results$Growth_Rate > 0,                       # If growth rate is positive
  log10(results$Growth_Rate),                    # Use log10 directly
  -log10(-results$Growth_Rate + 1e-6)           # For negative values, add a small constant to avoid log10(0) issues
)

# Calculate log_pvalue for use in the volcano plot
results$log_pvalue <- -log10(results$pvalue)

# Create a color_group variable to indicate significantly upregulated and downregulated genes
results$color_group <- ifelse(
  results$pvalue < 0.05,
  ifelse(
    results$Growth_Rate > 0,   # Upregulated gene
    "Red",                    # Upregulated
    "Blue"                    # Downregulated
  ),
  "Not Significant"            # Not significant
)

# Match gene names to results
results$GeneName <- gene_info$GeneName[match(results$Gene, gene_info$GeneID)]



# Assuming average_cpm has already been computed
exp_last_round <- exp_cpm_normalized[[length(exp_cpm_normalized)]]
exp_first_round <- exp_cpm_normalized[[1]]

#Optional: Background Substraction
control_last_round <- control_cpm_normalized[[length(control_cpm_normalized)]]
control_first_round <- control_cpm_normalized[[1]]

#Optional: Background Substraction
average_cpm <- rowMeans((exp_last_round - control_last_round) - (exp_first_round - control_first_round))

# Convert to log10 values and save to results
results$log10_CPM <- log10(average_cpm + 1e-6)  # Add a pseudocount to avoid log10(0)

results$abs_log10_CPM <- abs(results$log10_CPM)

# Filter highlight_genes to include only those with corresponding points in results
results <- results %>%
  filter(Gene %in% results$Gene & 
           !is.na(normalized_Growth_Rate) & 
           !is.na(abs_log10_CPM) & 
           !is.na(color_group) & 
           !is.na(log_pvalue))


results <- results %>%
  filter(abs_log10_CPM >= 3)  # remove rows that log10CPM is less than 3 which means CPM less than 1000

# Define highlight_genes based on color_group
highlight_genes <- results %>%
  filter(color_group %in% c("Red", "Blue"))


# Select top 30 genes by log10_CPM
top_genes_CPM <- highlight_genes %>%
  arrange(desc(log10_CPM)) %>%
  slice_head(n = 30)

# Select top 30 genes by -log10(p-value)
top_genes_pvalue <- highlight_genes %>%
  arrange(desc(log_pvalue)) %>%  # Sort by log_pvalue to get smallest p-values
  slice_head(n = 30)

# Combine the two sets of top genes
top_genes_combined <- bind_rows(top_genes_CPM, top_genes_pvalue) %>%
  distinct()  # Remove duplicates if any gene appears in both sets

# Adjust the legend settings and ensure plot uses space efficiently
volcano_plot <- ggplot(results, aes(x = normalized_Growth_Rate, y = log_pvalue)) +
  geom_point(aes(color = color_group, size = abs_log10_CPM), alpha = 0.6) +
  scale_color_manual(
    values = c("Red" = "red", "Blue" = "blue", "Not Significant" = "grey"),
    labels = c("Significantly Downregulated", "Not Significant", "Significantly Upregulated")
  ) +
  theme_minimal() +
  labs(
    x = "Log10 Growth Rate", y = "-log10(p-value)",
    color = "Gene Regulation", size = "Absolute Log10 CPM"
  ) +
  geom_text_repel(
    data = top_genes_combined,
    aes(label = ifelse(!is.na(GeneName) & GeneName != "", GeneName, Gene)),
    size = 3, box.padding = 0.3, point.padding = 0.5, arrow = arrow(length = unit(0.01, "npc")),
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
    plot.margin = margin(5, 5, 5, 5, "pt")  # Adjust margins if necessary
  ) 


# Output paths
plot_base_dir <- file.path(exp_base_path, "Routput")

# Ensure the directory exists or create it
if (!dir.exists(plot_base_dir)) {
  dir.create(plot_base_dir, recursive = TRUE)
}

# Set row names to values in the "Gene" column
rownames(differences_df_named) <- differences_df_named$Gene

# Optionally, remove the "Gene" column as it's now redundant
differences_df_named <- differences_df_named %>%
  select(-Variance) %>% select(-Gene)

# Save the modified data frame to a CSV file (add row names)
csv_file_path <- file.path(plot_base_dir, "differences_df_named_no_variance.csv")
write.csv(differences_df_named, csv_file_path, row.names = TRUE)



# After creating plots, ensure you save them into the existing directory
ggsave(file.path(plot_base_dir, "pca_plot.png"), plot = pca_plot, width = 8, height = 6)
ggsave(file.path(plot_base_dir, "line_plot.png"), plot = line_plot, width = 8, height = 6)
ggsave(file.path(plot_base_dir, "ttest_volcano_plot.png"), plot = volcano_plot, width = 12, height = 6)

