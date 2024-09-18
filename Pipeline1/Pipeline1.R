# List of required packages
required_packages <- c(
  "BiocManager", "Signac", "rtracklayer", "ggplot2", 
  "Rsubread", "DESeq2", "edgeR", "dplyr", "tidyr"
)

# Function to check if a package is installed and install it if not
install_missing_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% rownames(installed.packages())) {
        # If it's already installed just load it
        library(pkg, character.only = TRUE)
      } else {
        # If it's not installed via Bioconductor, use install.packages for CRAN packages
        if (pkg %in% c("ggplot2", "dplyr", "tidyr")) {
          install.packages(pkg)
        } else {
          # Use BiocManager for Bioconductor packages
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager")
          }
          BiocManager::install(pkg)
        }
      }
    }
  }
}

# Check and install missing packages
install_missing_packages(required_packages)


library(rtracklayer)
library(ggplot2)
library(Rsubread)
library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyr)

#Path to your newly generated Routput folder
plot_base_dir <- "/srv/scratch/z3546698/true/Routput"

# Read GTF file
gtf_file <- "/srv/scratch/z3546698/true/reference/Homo_sapiens.GRCh38.110.gtf"  # Please replace with the actual path
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

# Define paths
exp_base_path <- "/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119"
#exp_base_path <- "/srv/scratch/z3546698/true/Small_Molecule/JQ1/CoT/240302/"


control_base_path <- "/srv/scratch/z3546698/true/Small_Molecule/Biotin/T7MB-2/240421"
#control_base_path <- "/srv/scratch/z3546698/true/Small_Molecule/Biotin/CoT/240413"

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
    step4_dir <- file.path(base_dir, dir, "step3")
    
    # Construct the full file path to the expression counts file
    file_path <- file.path(step4_dir, paste0(dir, "_combined_expression_counts.txt"))
    
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

# Select the top 20 genes with the largest variance
top_genes <- differences_df_named %>%
  arrange(desc(Variance)) %>%
  head(20) %>%
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
  scale_color_manual(values = c("Control" = "red", "Experiment" = "blue")) +
  scale_fill_manual(values = c("Control" = "red", "Experiment" = "blue")) +  # Set fill color
  theme_classic() +
  labs(title = "Principal Component Analysis",
       x = "Principal Component 1",
       y = "Principal Component 2")

ggsave(file.path(plot_base_dir, "pca_plot.png"), plot = pca_plot, width = 8, height = 6)



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
  geom_line() +
  geom_point() +
  scale_color_manual(values = colors) +  # Use custom colors
  labs(
    title = "Top Genes with the Largest Change Across Timepoints (exp)",
    x = "Timepoints",
    y = "Normalized Gene Counts",
    color = "Gene"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right"  # Or "bottom", "none", adjust as needed
  )


ggsave(file.path(plot_base_dir, "line_plot.png"), plot = line_plot, width = 8, height = 6)



first_timepoint <- timepoints[1]

last_timepoint <- timepoints[length(timepoints)]

differences_df_named$Growth_LR_FR <- differences_df_named[last_timepoint] - differences_df_named[first_timepoint]

# Extract data for the first round & last round

exp_FR <- exp_cpm_normalized[[first_timepoint]]

exp_LR <- exp_cpm_normalized[[last_timepoint]]

# Ensure data are columns of a matrix
exp_FR <- as.matrix(exp_FR)
exp_LR <- as.matrix(exp_LR)

results <- data.frame(Gene=rownames(exp_LR), pvalue=NA)



for (i in 1:nrow(exp_FR)) {
  # Extract expression values for each gene in first round and last round
  fr_values <- exp_FR[i, ]
  lr_values <- exp_LR[i, ]
  
  # Perform unpaired t-test
  t_test_result <- t.test(fr_values, lr_values, paired = FALSE)
  
  # Extract and save the p-value
  results$pvalue[i] <- t_test_result$p.value
}


matching_genes <- intersect(rownames(differences_df_named), results$Gene)

# Extract matching p-values
matching_pvalues <- results$pvalue[results$Gene %in% matching_genes]

# Save matching p-values into differences_df_named
differences_df_named$pvalue[rownames(differences_df_named) %in% matching_genes] <- matching_pvalues


differences_df_named$Growth_LR_FR <- differences_df_named[last_timepoint] - differences_df_named[first_timepoint]

# Log transform the growth rate to manage outliers and improve distribution

library(ggplot2)
library(dplyr)
                  
epsilon <- 1e-6

# Applying signed log transformation
differences_df_named$Normalized_Growth_LR_FR <- ifelse(
  differences_df_named$Growth_LR_FR > 0,
  scale(log10(differences_df_named$Growth_LR_FR + 1)),
  scale(-log10(-differences_df_named$Growth_LR_FR + 1))
)


differences_df_named$log_pvalue <- -log10(differences_df_named$pvalue)


differences_df_named$color_group <- ifelse(
  differences_df_named$pvalue < 0.05,
  ifelse(
    differences_df_named[[last_timepoint]] > 2 * differences_df_named[[first_timepoint]],
    "Red",
    ifelse(
      differences_df_named[[last_timepoint]] < 0.5 * differences_df_named[[first_timepoint]], 
      "Blue", 
      "Not Significant"
    )
  ),
  "Not Significant"
)

highlight_genes <- differences_df_named %>%
  filter(color_group %in% c("Red", "Blue"))


# Create the color_group variable under the new conditions.
differences_df_named$color_group <- ifelse(
  differences_df_named$pvalue < 0.05,
  ifelse(
    differences_df_named[[last_timepoint]] > 2 * differences_df_named[[first_timepoint]],
    "Red",
    ifelse(
      differences_df_named[[last_timepoint]] < 0.5 * differences_df_named[[first_timepoint]], 
      "Blue", 
      "Not Significant"
    )
  ),
  "Not Significant"
)

highlight_genes <- differences_df_named %>%
  filter(color_group %in% c("Red", "Blue"))

# Define a scaling function using a power transformation
scale_size_adjust <- function(x) {
  scale_min <- 0.1  # Small minimum size for near-zero values
  scale_max <- 10   # Large maximum size for large values
  power <- 2.5      # Power factor to intensify scaling
  
  # Normalize x between 0 and 1
  norm_x <- (abs(x) - min(abs(x))) / (max(abs(x)) - min(abs(x)))
  
  # Use power transformation for scaling
  scaled_size <- scale_min + (scale_max - scale_min) * (norm_x ^ power)
  return(scaled_size)
}

# Apply the function and plot
volcano_plot <- ggplot(differences_df_named, aes(x = Normalized_Growth_LR_FR, y = log_pvalue)) +
  geom_point(aes(color = color_group, size = scale_size_adjust(Normalized_Growth_LR_FR)), alpha = 0.6) +
  scale_color_manual(values = c("Red" = "red", 
                                "Blue" = "blue",
                                "Not Significant" = "darkgrey"),
                     labels = c("Significantly Downregulated", 
                                "Not Significant", 
                                "Significantly Upregulated")) +
  theme_minimal() +
  labs(
    x = "Growth Rate",
    y = "-log10(p-value)",
    title = "Volcano Plot of Expression",
    color = "Gene Regulation",
    size = "Scaled Size"
  ) +
  geom_text(data = highlight_genes, aes(label = GeneName),
            size = 4, vjust = -0.5, hjust = 0.5, check_overlap = TRUE)

                  
ggsave(file.path(plot_base_dir, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)

