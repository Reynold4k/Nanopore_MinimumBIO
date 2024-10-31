

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


## Setup Instructions
   
## Step1 Download genome reference via:
Genome reference file is used for alignment from your sequence files to the reference.

Here's an illustration of downloading and decompressing hg38 human reference genome:

```bash
# Install the wget module
apt install wget

# Download the human reference genome file: hg38.fa
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

#Install the gunzip decompressor:
pip install gunzip

# Decompress the human reference genome file
gunzip hg38.fa.gz

```

## Step2 Index your reference fasta files, the successful implement of this step would ensure you won't meet errors in the later steps

It usually takes you a while.

```bash
#Download BWA software via:
apt install bwa

#Index using BWA
bwa index hg38.fa
```

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




### Script Usage Instruction

## Action 0 Installation of required modules

```bash
#Porechop
apt install porechop

#Nanoplot
pip install NanoPlot

#Samtools
apt install samtools

#Subread
apt install subread

#Seqkit
apt install seqkit

apt install nano
```


Set the input and reference data paths in the script:

## Action 1 Modify the preprocessing script(e.g. porechop_preprocessing.sh)
Do the same modification to the pod5 files processing scripts.

For fastq files: porechop_preprocessing.sh
For Pod5 files: pod5_preprocessing_with_genecounts.sh


- **Modify the `FOLDER` variable**: Set it to the directory containing your sequencing files.

If you aren't sure about your location, you can print your location just by running:
```bash
pwd
```

For windows users only, you may need to convert the script format if you meet errors like "command not found":

```bash

dos2unix porechop_preprocessing.sh

```


You may have different ways of opening the script:
1. Nano the_name_of_the_script
Please check this link to know more about Nano command: https://ioflood.com/blog/nano-linux-command/#:~:text=To%20use%20the%20nano%20command,will%20create%20it%20for%20you
  
2. Directly use text editor to open via system operation portal
Use Control + S to save your script aftering editing

```bash
# What you need to modify:
# 1.Change the line 19 FOLDER to your/pod5/path
# 2.Change the line 25 REFERENCE to your/reference/path
# 3.Change the line 27 ANNOTATION to your/ANNOTATION/path


FOLDER="the outout of your working directory having the same structure as below"

For example,
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders
```
![image](https://github.com/user-attachments/assets/c507eded-0724-438a-afc2-74a3484a5785)



## Action 2 Update gene reference and annotation file path in the pre processing script

- **Update the `REFERENCE` and `annotation` paths**: Point these to your reference genome and annotation files, respectively.

For example,

```bash
REFERENCE="/mnt/d/hg38/hg38.fa"
annotation="/mnt/d/hg38/hg38.ensGene.gtf"
```

Ensure your directory structure allows the script to find input files and output directories correctly as described in the script comments.

```bash
#If you want to accelerate your porechop that you supposed to have a huge sequences, then you need to add threads:

# change "-t 16" to the maximum number of CPU's you booked as this step can be slow
porechop -t 16 -i "$merged_file" -o "$output_file" > "$log_file" 2>&1

```

### Before running, if you're having another known in-frame check sequence different from: GATCCGAATTCN

Please replace the seqkit line in the porechop_preprocessing.sh with:

```bash

seqkit replace -p "^((?:.*?\n){3}).*?(your_in_frame_check_sequence)(\n.*)" -r '$1$2$3' -o "$seqprocessed_file" "$merged_file"

For example:
seqkit replace -p "^((?:.*?\n){3}).*?(GATCATTACTGAGCTATAGCTCATGCGGCCGC)(\n.*)" -r '$1$2$3' -o "$seqprocessed_file" "$merged_file"

```

## Action 3 Save the edited script and execute script

Use control + s in the keyboard to save the scripts

- **Execute the scripts for their corresponding data

For example, for nanopore fastq data, run command:


```bash 

bash ./porechop_preprocessing.sh

```

### Porechop Usage and Configurations

While trimming adapters from sequencing reads using Porechop, you may want to optimize performance and handle any potential errors effectively. Below are guidelines and options specific to running Porechop:

#### What If There Are Errors?

- **Log Files**: During execution, standard messages are redirected to a log file, while errors are captured separately for easier troubleshooting. Review `porechop_error.log` for any error messages to understand and resolve issues.


If run successfully, you'll see:
![image](https://github.com/user-attachments/assets/208572f3-bbff-4e30-acc2-107beacb8476)


### Output


Quality Control Outputs: Assumes NanoPlot produces a directory for quality control results and checks for this directory.
(For details, please refer to the Nanoplot github page for the explanation of different plots:https://github.com/wdecoster/NanoPlot)

BAM Files: Validates the existence of BAM files in step1 directories.
Sorted and Marked BAM Files: Checks for sorted and marked BAM files in step2 directories.
FeatureCounts Outputs: Confirms the presence of expression count files in step3 directories."




## Action 4 Analysis of generated gene counts matrix: Updating 3 Paths in R and Bash Scripts

## Introduction to R and RStudio

R is a powerful programming language and software environment used for statistical computing and graphics. It's popular among statisticians and data scientists for data analysis and visualization. RStudio is an integrated development environment (IDE) that enhances the R programming experience.

### How to Download and Install R and RStudio

### Step 1: Download and Install R

1. **Visit the CRAN Website:**
   - Go to the [CRAN R Project website](https://cran.r-project.org/).

2. **Choose Your Operating System:**
   - Click on the link for your operating system (Windows, macOS, or Linux).

3. **Download R:**
   - Follow the link to download the latest version of R for your OS.
   - **For Windows:** Click "Download R for Windows" and then "base" for the installer.
   - **For macOS:** Click on the "R-4.x.x.pkg" (version number may vary).

4. **Install R:**
   - Run the downloaded installer and follow the installation steps.

### Step 2: Download and Install RStudio

1. **Visit the RStudio Website:**
   - Go to the [RStudio Download page](https://www.rstudio.com/products/rstudio/download/).

2. **Choose RStudio Desktop:**
   - Click the "Download" button for RStudio Desktop (Open Source License).

3. **Download RStudio:**
   - Choose the installer for your operating system and download it.

4. **Install RStudio:**
   - Run the installer and follow the instructions to install RStudio.

### Step 3: Getting Started with R and RStudio

1. **Open RStudio:**
   - Launch RStudio from your applications folder or start menu.

2. **Familiarize Yourself with the Interface:**
   - RStudio's interface is divided into multiple panes: console, script editor, environment/history, and files/plots/packages/help. Explore these to understand their functions.

3. **Write Your First R Script:**
   - In the script editor, type: `print("Hello, World!")`.
   - Run the script using the "Run" button or press `Ctrl + Enter` (Cmd + Enter on macOS).

4. **Explore R Packages:**
   - Extend R's functionality by installing packages with: `install.packages("package_name")`.

5. **Learn R Basics:**
   - Get started with R's basic data structures like vectors, matrices, data frames, and lists, and learn basic operations such as loops and functions.

## Additional Resources

- **R Documentation:** The official [R Documentation](https://www.rdocumentation.org/) is comprehensive and detailed.
- **Online Courses:** Platforms like Coursera, edX, and DataCamp offer structured R courses.
- **Communities and Forums:** Engage with communities on Stack Overflow or Reddit for help and support.


## After you downloaded Rstudio and installed R:

Step A. Change the to your exp and control folders
In your R scripts, make sure that `exp_base_path` and `control_base_path` paths match the formats of `FOLDER` path specified in the `porechop_preprocessing.sh` script. Below is an example of how you can configure these: 

```r
# Example Path configuration in the R script
exp_base_path <- "/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"
control_base_path <- "/mnt/d/Small_Molecule/JQ1/T7MB-2/240421"

```

Step B. Change the  ANNOTATION to your/ANNOTATION/path

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

Step C. Packages installation:

Open the linux portal and type:

```bash
# To install R
sudo apt install r-base

# To open R
R
```

Inside the R terminal you opened when typing R, run the following installation:

```r
# Updating Bioconductor and all necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

install.packages(c("Matrix", "MASS", "mgcv", "ggplot2", "dplyr", "tidyr"))

# Bioconductor packages
BiocManager::install(c("DelayedArray", "SummarizedExperiment", "DESeq2", "rtracklayer"))

```

Update all/some/none? [a/s/n]: a

### !!! Answer a to make sure all your packages are up-to-date

If you never installed packages in your PC, it would take a while to install everything, be patient.

Then quit R via, enter yes to save the workspace:

```r


```

### After you modified the path directory in the Analysis.R, run the R script in the linux portal(if windows users using Ubuntu, mac users using terminal, please refer to the guideline in the main menu) through:
   
   ```bash
    Rscript Analysis.R
   ``` 
### If you prefer to manually run R scripts using RStudio instead of the command line, follow these steps:

1. **Open RStudio:**

   - Launch RStudio from your applications folder or start menu.

2. **Create a New Script File:**

   - Click on `File` > `New File` > `R Script` to open a new script editor.

3. **Copy and Paste the Code:**

   - Copy your R script content (from `Analysis.R` or another file) and paste it into this new script editor.

4. **Save the Script:**

   - Save the script by selecting `File` > `Save` or pressing `Ctrl + S` (`Cmd + S` on macOS). Name the file appropriately, such as `MyAnalysis.R`.

5. **Run the Script:**

   - **Step A: Check Environment and Dependencies:**

     - Ensure all necessary packages are installed. Install any missing packages using:
       ```r
       install.packages(c("Matrix", "MASS", "mgcv", "ggplot2", "dplyr", "tidyr"))

       # Bioconductor packages
       BiocManager::install(c("DelayedArray", "SummarizedExperiment", "DESeq2", "rtracklayer"))
       ```

   - **Step B: Execute Code Lines or Chunks:**

     - To run specific lines or sections of the script, highlight the desired code portion and press `Ctrl + Enter` (`Cmd + Enter` on macOS), or click the "Run" button in the script editor toolbar.

   - **Step C: Run the Entire Script:**

     - To execute the entire script at once, click `Source` at the top-right of the script editor, or use the shortcut `Ctrl + Shift + S` (`Cmd + Shift + S` on macOS).

6. **View Results:**

   - Outputs from your script are shown in the console pane at the bottom of RStudio. Review this pane for any messages or errors and adjust your code as necessary.



# R Script Adjustment Guide


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









