

```bash
 " __  __ ___ _   _ ___ __  __ _   _ __  __   _     _       "
 "|  \/  |_ _| \ | |_ _|  \/  | | | |  \/  | | |__ (_) ___  "
 "| |\/| || ||  \| || || |\/| | | | | |\/| | | '_ \| |/ _ \ "
 "| |  | || || |\  || || |  | | |_| | |  | | | |_) | | (_) |"
 "|_|  |_|___|_| \_|___|_|  |_|\___/|_|  |_| |_.__/|_|\___/ "


```

# Pipeline1

## Features
1. **In-Frame Check**: Utilizes Seqkit for only keeping the sequence in frame and removing all other (backwards/out of frame)sequences.
2. **Trimming**: Utilizes Porechop for removing adapters from sequencing reads.
3. **Quality Control**: Leverages NanoPlot for assessing read quality.
4. **Alignment**: Aligns reads to a reference genome using BWA.
5. **Sorting and Deduplication**: Sorts and marks duplicate reads using Samtools.
6. **Gene Feature Counting**: Uses featureCounts for quantifying gene expression.

## Prerequisites

Ensure your system meets the following requirements before running the pipeline:

- **Operating System**: Linux or compatible environment.

# For Windows Users Ubuntu from microsoft store is highly recommended:
- You could also obtain this software via: https://ubuntu.com/

![image](https://github.com/user-attachments/assets/9b400296-fc30-4dbd-929c-7c9bdafaf438)

## Guide for Windows Users Transitioning to Ubuntu

## Understanding File Paths
- **Path Delimiters:**  
  Unlike Windows, which uses backslashes (`\`) for file paths, Ubuntu utilizes forward slashes (`/`).  
  Example: `C:\Users\YourName\Documents` becomes `/home/YourName/Documents` in Ubuntu.

- **Case Sensitivity:**  
  File and directory names in Ubuntu are case-sensitive.  
  `Documents`, `documents`, and `DOCUMENTS` are considered different directories.

## Home Directory
- On Ubuntu, your personal files are stored in the home directory, typically located at `/home/YourName`. 
  This is similar to `C:\Users\YourName` on Windows.

## Accessing Drives
- Windows drives (C:, D:, etc.) are mounted in the `/mnt` or `/media` directory in Ubuntu.  
  For example, your C: drive might be accessible under `/mnt/c`.

## Hidden Files and Directories
- Files and folders prefixed with a dot (`.`) in Ubuntu are hidden by default.  
  To view them in the file manager, press `Ctrl + H`.

## Permissions
- Ubuntu enforces file permissions more strictly than Windows.  
  You might need to modify permissions using commands like `chmod` or change the owner with `chown` for certain tasks.

## Using the Terminal
- The terminal is a powerful tool in Ubuntu, used for a variety of tasks.  
  Get familiar with basic commands like `ls` (list), `cd` (change directory), and `cp` (copy).

## Installing Software
- Unlike Windows, software in Ubuntu is often installed via package managers like `apt` (Advanced Package Tool).  
  You might also use software repositories like the Ubuntu Software Center.

## File Extensions
- Ubuntu does not rely on file extensions to identify file types as strictly as Windows.  
  It often determines the file type by its content.

## Backup Your Data
- Before making any major changes, always ensure your data is backed up.  
  Tools like `rsync` can be invaluable for maintaining backups in Ubuntu.


# For Mac Users Click the Launchpad icon in the Dock, type Terminal in the search field, then click Terminal.
# In the Finder , open the /Applications/Utilities folder, then double-click Terminal.

    
## Beginner's Guide to Using Terminal and Command Line on macOS

## Introduction to Terminal
- **Accessing Terminal:**  
  Terminal is a built-in application in macOS, found under `Applications > Utilities > Terminal`.

- **Basic Interface:**  
  Terminal provides a command line interface where you can type commands to perform various tasks.

## Navigating the File System
- **Current Directory:**
  - Use `pwd` to print the current working directory.

- **Listing Files and Directories:**  
  - Use `ls` to list files and directories.  
  - Use `ls -la` to include hidden files and detailed information.

- **Changing Directory:**  
  - Use `cd [directory_name]` to navigate to a different directory.  
  - Use `cd ..` to go up one directory level.

## Managing Files and Directories
- **Creating Directories:**  
  - Use `mkdir [directory_name]` to create a new directory.

- **Creating Files:**  
  - Use `touch [file_name]` to create a new, empty file.

- **Copying Files:**  
  - Use `cp [source] [destination]` to copy files or directories.

- **Moving/Renaming Files:**  
  - Use `mv [source] [destination]` to move or rename files.

- **Deleting Files and Directories:**  
  - Use `rm [file_name]` to delete files.  
  - Use `rm -r [directory_name]` to delete directories and their contents.

## Editing Files
- **Using Nano Editor:**  
  - Use `nano [file_name]` to edit files directly in Terminal.

## Permissions
- **Changing Permissions:**  
  - Use `chmod [permissions] [file_name]` to change file or directory permissions.

- **Changing Ownership:**  
  - Use `chown [user] [file_name]` to change file ownership.

## Searching and Finding Files
- **Search with grep:**  
  - Use `grep [search_term] [file_name]` to search for a term within a file.

- **Finding Files:**  
  - Use `find [directory] -name [file_name]` to search for files by name.

## Useful Tips
- **Auto-Completion:**  
  - Use the `Tab` key to auto-complete file and directory names.

- **Command History:**  
  - Use the `Up` and `Down` arrow keys to cycle through command history.

- **Canceling Commands:**  
  - Use `Ctrl + C` to cancel an ongoing command or process.


**Reference and Annotation Files**:
- Reference genome in FASTA format (e.g., `hg38.fa`).
- Gene annotation file in GTF format (e.g., `hg38.ensGene.gtf`).


## Setup Instructions

## Sometimes the rclone mount would not be able to mount all the data correctly, just make sure that rclone mount is working properly.

   
## Download hg38 reference and annotation files via:

```bash

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

## Index your reference fasta files, this step must be conducted on your current portal to continue the next step
If you want to use a different reference genome for the `bwa index` command instead of `hg38.fa`, you can follow these steps:

1. **Download or Prepare Your Reference Genome**: First, ensure that you have the FASTA file for the reference genome you'd like to use. This could be any genome such as `mm10.fa` (mouse), `GRCh37.fa` (older human reference), or even a custom genome assembly. The file should be in FASTA format, which usually has a `.fa` or `.fasta` extension.

2. **Modify the Command**: Replace `hg38.fa` with the path to your reference genome. For example, if you want to index the `GRCh37.fa` reference genome, the command would look like this:

```bash

bwa index /path/to/your/reference_genome.fa

#For example:
bwa index GRCh37.fa

```

This section of commands performs two main actions: downloading and decompressing a GTF (General Feature Format) file that contains gene annotation data from Ensembl. Here’s a detailed breakdown:

**Download the GTF File**:

```bash

wget http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz

```

**Decompress the GTF File**:


```bash
gunzip Homo_sapiens.GRCh38.110.gtf.gz
```
### Note: You can download GTF files for other species from the Ensembl FTP server or UCSC Genome Browser.
### Below are examples of how to find the corresponding paths:

 1. Ensembl FTP Server:
Visit the Ensembl FTP site: ftp://ftp.ensembl.org/pub/release-110/gtf/
To download a GTF for Mus musculus (mouse), use:
ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/
### Example path for Mus musculus: "/path/to/your/directory/Mus_musculus.GRCm39.110.gtf.gz"

2. UCSC Genome Browser:
Visit the UCSC Genome Browser downloads section: http://hgdownload.soe.ucsc.edu/downloads.html
To find a GTF for Danio rerio (zebrafish), browse to:
http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/genes/
Example path for Danio rerio: "/path/to/your/directory/danRer11.refGene.gtf.gz"

Make sure to decompress (.gz) files if needed and adjust the path accordingly.




### Usage Instruction

Set the input and reference data paths in the script:

## Step1 Modify the preprocessing script(e.g. porechop_preprocessing.sh)
Do the same modification to the pod5 files processing scripts

For fastq files: porechop_preprocessing.pbs
For Pod5 files: pod5_preprocessing_with_genecounts.pbs


- **Modify the `FOLDER` variable**: Set it to the directory containing your sequencing files.

If you aren't sure about your location, you can print your location just by running:
```bash
pwd
```

```bash
# What you need to modify:
# 1.Change the line 25 FOLDER to your/pod5/path
# 2.Change the line 26 REFERENCE to your/reference/path
# 3.Change the line 27 ANNOTATION to your/ANNOTATION/path


FOLDER="the outout of your working directory having the same structure as below"

For example,
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders
```
![image](https://github.com/user-attachments/assets/c507eded-0724-438a-afc2-74a3484a5785)



## Step2 Update gene reference and annotation file path in the pre processing script

- **Update the `REFERENCE` and `annotation` paths**: Point these to your reference genome and annotation files, respectively.

For example,
REFERENCE="/mnt/d/hg38/hg38.fa"
annotation="/mnt/d/hg38/hg38.ensGene.gtf"

If you are using UNSW KATANA, the working directory format may differ

Ensure your directory structure allows the script to find input files and output directories correctly as described in the script comments.



## Step3 Save the edited script and execute script

- **Execute the scripts for their corresponding data

For example, for nanopore fastq data, run command:

```bash
#If you want to accelerate your porechop that you supposed to have a huge sequences, then you need to add threads:

# change "-t 16" to the maximum number of CPU's you booked as this step can be slow
porechop -t 16 -i "$merged_file" -o "$output_file" > "$log_file" 2>&1

```
### Porechop Usage and Configurations

While trimming adapters from sequencing reads using Porechop, you may want to optimize performance and handle any potential errors effectively. Below are guidelines and options specific to running Porechop:

#### What If There Are Errors?

- **Log Files**: During execution, standard messages are redirected to a log file, while errors are captured separately for easier troubleshooting. Review `porechop_error.log` for any error messages to understand and resolve issues.

#### Configuration Options:

- **Log and Error Files**: The output and any errors will be logged as follows:
  ```bash
  # Standard output and errors
  --verbosity 1>> "$log_file" 2>> porechop_error.log


### Before running, if you're having another known in-frame check sequence different from: GATCCGAATTCN

Please replace the seqkit line in the porechop_preprocessing.sh with:

```bash

seqkit replace -p "^((?:.*?\n){3}).*?(your_in_frame_check_sequence)(\n.*)" -r '$1$2$3' -o "$seqprocessed_file" "$merged_file"

For example:
seqkit replace -p "^((?:.*?\n){3}).*?(GATCATTACTGAGCTATAGCTCATGCGGCCGC)(\n.*)" -r '$1$2$3' -o "$seqprocessed_file" "$merged_file"

```

```bash 

bash ./porechop_preprocessing.sh

```

If run successfully, you'll see:
![image](https://github.com/user-attachments/assets/208572f3-bbff-4e30-acc2-107beacb8476)




#### Functionality of the Script:

1. **Define the Main Directory**:
   - Set the `FOLDER` variable to specify the path of the main directory where your experimental folders are located.
   - Replace `"/path/to/your/folder/you/want/to/clear"` with the path to your working directory.

2. **Check for Expected Output Files**:
   - The script verifies the existence of specific output files generated during various processing steps:
     - **Trimmed Files**: Checks for the presence of trimmed FASTQ files in `step1` directories.
     - **Quality Control Outputs**: Assumes NanoPlot produces a directory for quality control results and checks for this directory.
     - **BAM Files**: Validates the existence of BAM files in `step2` directories.
     - **Sorted and Marked BAM Files**: Checks for sorted and marked BAM files in `step3` directories.
     - **FeatureCounts Outputs**: Confirms the presence of expression count files in `step4` directories.

3. **Flag for Completeness**:
   - The variable `all_outputs_exist` starts as `true`. If any expected output file is missing, it switches to `false`, and the script outputs the missing file(s).

4. **Conditional Deletion**:
   - If all outputs are confirmed to exist, the script deletes all contents within the `step1` to `step4` directories, thereby clearing intermediate files and folders.
   - If any expected output is missing, the script informs you and skips the deletion process to ensure no important data is lost.


## Output

The pipeline produces the following outputs:

- **Trimmed Fastq Files**: Located in the same directory as input with `_trimmed.fastq.gz` suffix.
- **Quality Reports**: Generated by NanoPlot in directories named as `<sample>_nanoplot`.
- **BAM Files**: Alignments stored as `<sample>.bam` files.
- **Sorted and Deduplicated BAM Files**: Named `<sample>_sorted_marked.bam`.
- **Expression Counts**: Gene expression data saved as `<folder_name>_expression_counts.txt`.



## step 4 Analysis of generated gene counts matrix: Updating 3 Paths in R and Bash Scripts

In your R script, make sure that `exp_base_path` and `control_base_path` paths match the `FOLDER` path specified in the `porechop_preprocessing.sh` script. Below is an example of how you can configure these: 

```r
# Example Path configuration in the R script
exp_base_path <- "/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"
control_base_path <- "/mnt/d/Small_Molecule/JQ1/T7MB-2/240421"

```

## Downloading GTF File on Linux

To download the `Homo_sapiens.GRCh38.112.gtf` file on a Linux system, use the following command in the terminal:

```bash
# If you haven't downloaded the gtf annotation file before, download it via:

wget ftp://ftp.ensembl.org/pub/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz -P /path/to/your/directory

```

## Then edit R scripts in the gtf path:

### This is included in the R scripts:

1.Change the  to your exp and control folders
 
```r
# The Path here should be in the same level as above Porechop_processing.pbs
exp_base_path <- "/srv/scratch/z3546698/tutorial/Small_Molecule/JQ1/Co/240217"
control_base_path <- "/srv/scratch/z3546698/tutorial/Small_Molecule/Biotin/Co/240217"

```

2.Change the  ANNOTATION to your/ANNOTATION/path
```r
gtf_file <- "/path/to/your/directory/Homo_sapiens.GRCh38.112.gtf.gz"
# Note: You can download GTF files for other species from the Ensembl FTP server or UCSC Genome Browser.
# Below are examples of how to find the corresponding paths:

# 1. Ensembl FTP Server:
# Visit the Ensembl FTP site: ftp://ftp.ensembl.org/pub/release-110/gtf/
# To download a GTF for Mus musculus (mouse), use:
# ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/
# Example path for Mus musculus: "/path/to/your/directory/Mus_musculus.GRCm39.110.gtf.gz"

# 2. UCSC Genome Browser:
# Visit the UCSC Genome Browser downloads section: http://hgdownload.soe.ucsc.edu/downloads.html
```


Then run through the whole R scripts and check the result plots.

```bash

#load r to a stable version has been passed the test
module load r/4.4.0
#To start R in the katana terminal
R
```

```r
# Updating Bioconductor and all necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

install.packages(c("Matrix", "MASS", "mgcv", "ggplot2", "dplyr", "tidyr"))

# Bioconductor packages
BiocManager::install(c("DelayedArray", "SummarizedExperiment", "DESeq2", "rtracklayer"))

```

### Answered "yes" to save workspace.
### Then quit r q()

### After you modified the path directory in the Analysis.R, run the R script through:
   
   ```bash
    Rscript Analysis.R
   ``` 

# R Script Adjustment Guide

This README provides instructions on adjusting the R script to modify PCA visualizations, color thresholds, and plot parameters for your RNA-seq data analysis. Follow these guidelines to customize your visual aspects effectively.

## Adjusting PCA Plot Parameters

To modify the appearance of the PCA plot, you can change various aesthetics and parameters.

- **Point Size in PCA Plot**: Change the size of the points to improve visibility in the PCA visualization.
  ```r
  geom_point(size = 5) +  # Adjust the point size
  ```

- **Color and Shape Aesthetics**: You can define colors and shapes based on sample grouping or timepoints.
  ```r
  aes(color = Group, shape = Timepoint)  # Change aes mappings as needed
  ```

- **Theme and Labels**: Customize the overall look by modifying the theme and axis titles for better clarity.
  ```r
  labs(title = paste(exp_name),
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  )
  ```

## Adjusting Color Thresholds

To change the criteria for color coding based on differential expression, edit the `case_when()` function in the script. This function determines which genes are considered significantly upregulated or downregulated.

### Color Assignment

- **Threshold for Color Coding**: Adjust the following code to set different cutoffs for coloring genes based on their expression changes.
  ```r
  mutate(
    color = case_when(
      log2foldchange < -1 ~ "blue",    # Modify this threshold for down-regulation
      log2foldchange > 1 ~ "red",      # Modify this threshold for up-regulation
      TRUE ~ "black"                   # Default color for no significant change
    )
  )
  ```

## Adjusting Plot Parameters

Adjusting plot aesthetics involves setting graphical parameters in ggplot functions. Here’s how you can change various aspects:

### Line Plot Parameters

- **Line and Point Size**: Change the `size` parameter in `geom_line()` and `geom_point()` to control how thick lines and how big the points appear in your line plots.
  ```r
  geom_line(size = 1.2)  # Increase or decrease for thicker or thinner lines
  geom_point(size = 3)   # Adjust for larger or smaller points
  ```

- **Text and Aesthetic Themes**: Utilize the `theme()` function to modify text sizes and positioning according to your presentation needs.
  ```r
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1),  # Adjust 'size' and 'angle'
    plot.title = element_text(size = 14, face = "bold"),           # Modify 'size' for the title
    axis.title = element_text(size = 20),                          # Change 'size' for axis titles
    legend.title = element_text(size = 18),                        # Update 'size' for legend titles
    legend.text = element_text(size = 14)                          # Set 'size' for legend text
  )
  ```

### Volcano Plot Parameters

- **Point Size and Alpha (Opacity)**: These are set in `geom_point()` where you can control the size and transparency of the points in the volcano plot.
  ```r
  geom_point(aes(color = color), alpha = 0.5, size = 5)  # Adjust 'size' and 'alpha'
  ```

- **Y-axis Limits and Theme**: Set y-axis limits using `ylim()` and adjust thematic elements using `theme_minimal()` and `theme()`.
  ```r
  ylim(1, 6)  # Modify y-axis bounds if necessary
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
  ```


If running successfully, you'll find the plot in the path of "exp_base_path" and see:

![image](https://github.com/user-attachments/assets/55010a3a-9f57-4faf-9632-379250f51131)









