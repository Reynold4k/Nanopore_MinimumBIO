
# Bioinformatics Data Processing Pipeline and Gene Coverage Analysis

## Overview of the Bioinformatics Data Processing Pipeline

This pipeline script is designed to process high-throughput sequencing data by executing a series of automated data processing steps on specified folders containing FASTQ files. The main functionalities include:

1. **Merging and Trimming FASTQ Files**: It merges all FASTQ files within each sample directory into a single file and trims the merged file using Porechop to remove low-quality sequences and adapter sequences.

2. **Quality Control**: Utilizes NanoPlot to assess the quality of the trimmed sequences, ensuring the reliability of the data.

3. **In frame check**: Utilizes Blastx to assess coding potentials of the data, removing sequences that out of frame like backwards.

4. **Alignment and Coverage Analysis**: Leverages Minimap2 to align the trimmed sequences to a reference genome, generating BAM files and calculating read coverage, which is then output to `coverage.txt`.

5. **Position Extraction and BED File Generation**: Extracts relevant gene positions from the coverage data and generates a BED file for easier visualization and downstream analysis.

6. **Sequence Extraction**: Based on the generated BED file, extracts specified regions from the reference sequences and saves them in FASTA format for further analysis and applications.

These scripts provides an automated workflow for bioinformatics research, significantly increasing data processing efficiency while ensuring accuracy and consistency in the results.

### Overall of all the scripts you would go through here is shown in order below:

pipeline1.5.sh -> Compare.sh -> Analysis.R

## Part1
### Step 0: Installation of several modules you would need to use in the later steps

If you've went through the pipeline1, then you only need to run:
```bash
#Install minimap2
apt install minimap2
```

or otherwisely:
```bash

# Install bc
apt install bc

# Install Python, if you didn't:
sudo apt install python3

# Install the wget module
apt install wget

#Install the gunzip decompressor:
pip install gunzip

#bwa
apt install bwa

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

In-frame check requires a protein database, here're the preparation steps:

### Downloading from UniProt

UniProt provides comprehensive sequence datasets of proteins. Here’s how to download and set up the human protein database for BLAST:

#### Get the Dataset

Use `wget` to fetch the FASTA file of all proteins from UniProt:

```bash
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```

### Filter for Human Proteins

After downloading, filter for human entries if the dataset contains multiple species. This can be done using:

```bash
zcat uniprot_sprot.fasta.gz | grep -A 1 '^>.*OS=Homo sapiens' > human_proteins.fasta
```

### Unzip and Prepare for BLAST
Unzip the compressed file and format it for use with BLAST:
```bash
gunzip uniprot_sprot.fasta.gz

#Install blast module on your pc for blast usage
apt-get install ncbi-blast+

makeblastdb -in human_proteins.fasta -dbtype prot -out human_proteome_db

```

### Index your reference file

You need to index your reference file before you're running everything.

```bash
#For example, for T7-pep:
bwa index Updated_T7_Pep_Ref_fixed.fasta
```

### Step 1: Modifying and Executing `pipeline1.5.sh`

In conclusion, what is going on in the pipeline1.5.sh is to go through every rounds in your specified experimental folder and load all the fastq first, then it would be merging all the fastq files into one to accelerate the analysing process, more details are available in the start of the script.

The software I recommend you to use is VS Code, you could download it via:https://code.visualstudio.com/.

So let's start!

To run the `pipeline1.5.sh` script, follow these steps:

1. **Open the Script**: Open the file `pipeline1.5.sh` in a text editor, or in commandline, use:
   
```bash
nano pipeline1.5.sh
```
Here's a complete guide for using nano command in the commandline:
https://ioflood.com/blog/nano-linux-command/#:~:text=To%20use%20the%20nano%20command,will%20create%20it%20for%20you.

2. **Modify Paths**:
    - Change FASTQ_FOLDER: Locate the line that defines `FASTQ_FOLDER` and replace it with the path to your folder containing the FASTQ files. For example:
      ```bash
      #You can add as many folder as you want as long as they all have the same reference file:
        FASTQ_FOLDER=(
            "/mnt/d/Bait_Glue/CRBN/glue/TON/240427"
            "/mnt/d/Bait_Glue/CRBN/MB014/TON/230827")
       ```
      This will be the folder path like below:
      
      ![image](https://github.com/user-attachments/assets/9cfb1d87-747c-4e3d-8f47-ec78c661b40b)

    - Change REFERENCE: Locate the line that defines `REFERENCE` and update the path to your reference file. For example:
      ```bash
      REFERENCE="/path/to/your/reference/file"
      ```
      
    - Change PROTEIN_DB: Locate the line that defines `PROTEIN_DB` and update the path to your protein database file, you just created it using "makeblastdb" command. For example:
      ```bash
      PROTEIN_DB="/mnt/d/human_proteome_db"
      ```      

![image](https://github.com/user-attachments/assets/942f7914-ddd5-4f75-b029-0f80e461b9cf)


3. **Set Up the Environment**: At the beginning of the script, add commands to load the required modules and set up a Conda environment. Follow these steps:


    
    - **Download miniconda**: 
    ```bash
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh
    source ~/.bashrc
    conda init
    ```

    - **Create the Conda Environment**: Create a Conda environment in your home directory.
    ```bash
    conda create -n my_conda python=3.10
    ```
    
    If you're unsure about your current directory path, using the below command:

    ```bash
    pwd
    ```
    
    Use the environment name from the directory where you want it:

    ```bash
    conda create --prefix my_conda python=3.10
    ```

    - **Activate the Conda Environment**: Activate the created Conda environment.
    ```bash
    conda activate my_conda
    ```
    Or, if you used a specific path to create the environment, you can activate it with:
    ```bash
    conda activate /your/current/directory/my_conda
    ```
    
    - **Install Required Packages**: Install the necessary Python packages with pip.
    ```bash
    pip3 install pandas matplotlib

    #Check if the installation is successful:
    which pandas
    which matplotlib
    ```
4. **Convert format**: If you're Windows users, you may need to convert the format of the scripts before running it, via:

   ```bash
   dos2unix pipeline1.5.sh
   ```

5. **Execute the Script**: Once you have modified the paths and set up the environment, save the file and submit it to your job scheduler. For example, use the following command:

    ```bash
    bash ./pipeline1.5.sh
    ```



### Expected Output

The output of the `pipeline1.5.sh` script will include:
- A merged and trimmed FASTQ file.
- A Blastx folder with in-frame check information
- Quality assessment report from NanoPlot.
- BAM files generated through the alignment process.
- Coverage data saved in `coverage.txt`.
- A BED file containing relevant gene positions.
- FASTA format files for specified regions.

These outputs will be stored in the specified output directories.

---


## Part2 Overview of the Gene Coverage Analysis Script `Compare.sh`

The main purpose of the `Compare.sh` script is to analyze the gene coverage between experimental and control groups, calculate differential coverage, and generate a corresponding horizontal bar plot for visualization. The specific steps are as follows:

1. **Load Required Modules**: Use the `module load` command to load the necessary Python and computational tools.

2. **Define Paths**: Specify the paths for the experimental group, control group, and reference files.

3. **Find the Latest Experimental Directory**: Retrieve the latest experimental and control round directories based on the folder structure.

4. **Check Coverage Files**: Confirm that the coverage files exist.

5. **Read Coverage Data**: Read the coverage data for the control and experimental groups into associative arrays.

6. **Calculate Differences**: Compare the coverage between the control and experimental groups, calculate the difference, and write the results to an output file.

7. **Generate Visualization Charts**: Use Python to create a horizontal bar plot that displays the coverage for the experimental group.

### Modifying and Executing `Compare.sh`

Function introduction: As you may find out yourself, if you don't really have an annotation gtf file rather than your personalized reference fasta, what you could do is to do the similar calculation on those data you could have at this stage. So Compare.sh compare the difference of the coverage files you generated during step2 between experimental and control groups when you specify them into the correct directory.

So let's continue!

To run the `Compare.sh` script, follow these steps:

1. **Open the Script**: Open the file `Compare.sh` in a text editor, or whatever you can to edit the text, if you're a person of super natural, please use the command:

Before running the script, you would be required to fully understand the functions of Nano command through, including how to edit and save it:
https://ioflood.com/blog/nano-linux-command/#:~:text=To%20use%20the%20nano%20command,will%20create%20it%20for%20you.

   ```bash
    Nano Compare.sh
   ```

2. **Modify Paths**:
    - Change EXPERIMENTAL_FOLDER: Locate the line that defines `EXPERIMENTAL_FOLDER` and replace it with the path to your experimental group's FASTQ parent folder. For example:
      ```bash
      EXPERIMENTAL_FOLDER="/path/to/your/exp/fastq"
      ```
    - Change CONTROL_FOLDER: Locate the line that defines `CONTROL_FOLDER` and replace it with the path to your control group's FASTQ parent folder. For example:
      ```bash
      CONTROL_FOLDER="/path/to/your/control/fastq"
      ```

![image](https://github.com/user-attachments/assets/5a2f7537-9bba-4d75-9802-6b22a9ee7991)


3. **Execute the Script**: Once you have modified the paths, save the file and submit it to your job scheduler. For example, use the following command:

    ```bash
    bash ./compare.sh
    ```


   
### Part3 Analysis: Once you have executed the scripts above successfully, analyse it using the R script (Analysis.R) following command:

## If you have a personalized library with gene names already existed in the reference file, please use "Analysis_withoutidmapping.R"
   
## What you what to modify:

## 1.Change the EXPERIMENTAL_FOLDER to your/exp/fastq parent path

```r
# Set the experimental folder path
EXPERIMENTAL_FOLDER <- "/mnt/d/Bait_Glue/VHL/MB012/TON/230827"
```

## 2.Change the id_mapping to your/id_mapping_file, you could save the id_mapping file in your previous reference folder if you want to:

Please note, "Analysis_withoutidmapping.R" does not require the id_mapping object below:

```r
id_mapping <- read.table("path/to/your/reference/idmapping_2024_10_01.tsv", 
                         header = TRUE, 
                         sep = "\t", 
                         stringsAsFactors = FALSE, 
                         fill = TRUE, 
                         quote = "",  
                         comment.char = "") 
```


## 3.Install the R packages that are not available:

Open R shell or Rstudio, you can find guideline here:
https://education.rstudio.com/learn/beginner/

```r
#install packages in the R portal via commands:

install.packages('gglpot2')
install.packages('dplyr')
install.packages('readr')
```
### If you're using R shell(similar to commanline), remeber that:

You open R in commanline via the command:
```r
R
```

Answered "yes" to save workspace.
Then quit r q()

### After you modified the path directory in the Analysis.R, run the R script directly in the R studio or through:
   
   ```bash
    Rscript Analysis.R
   ```

# R Script Adjustment Guide

This README provides instructions on adjusting the R script to modify plot appearances and color determination criteria for your differential coverage data analysis. It guides you on which parts of the code to edit to achieve desired changes.

## Adjusting Plot Appearance

To modify the visual aspects of your plots, you can tweak various aesthetic parameters in your ggplot2 calls:

### Line Plot Appearance

- **Line and Point Size**: Adjust the `geom_line()` and `geom_point()` functions to change the thickness of the lines and size of the points.
  ```r
  geom_line(size = 1.2)  # Modify 'size' for line thickness
  geom_point(size = 3)   # Modify 'size' for point size
  ```

- **Text and Theme**: You can alter font sizes and styles using the `theme()` function. This includes axis text, titles, and legend configuration.
  ```r
  theme(
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1),  # Adjust 'size' and 'angle' for x-axis text
    plot.title = element_text(size = 14, face = "bold"),           # Adjust 'size' for plot titles
    axis.title = element_text(size = 20),                          # Adjust 'size' for axis titles
    legend.title = element_text(size = 18),                        # Adjust 'size' for legend titles
    legend.text = element_text(size = 14)                          # Adjust 'size' for legend text
  )
  ```

### Volcano Plot Appearance

- **Filtering threshold of log10CPM**: change the number to your interested threshold, 3 means 10^3 = 1000.
```r
df <- df %>%
  filter(log10CPM >= 3)
```

- **Point Size and Transparency**: To change point size and transparency on the volcano plot, modify the appropriate parameters in `geom_point()`.
  ```r
  geom_point(aes(color = color), alpha = 0.5, size = 5)  # Adjust 'size' for point size and 'alpha' for transparency
  ```

- **Theme and Axis Limits**: For consistent styling, use `theme()` to manage text sizes and style. Use `ylim()` to adjust y-axis bounds.
  ```r
  ylim(1, 6)  # Adjust y-axis limits as needed
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 24),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
  ```

## Adjusting Color Determination

- **Color Assignment Based on Fold Change**: The code uses `case_when()` to determine colors based on log2 fold change. Adjust thresholds to meet different significance levels.
  ```r
  mutate(
    color = case_when(
      log2foldchange < -1 ~ "blue",    # Threshold for 'blue'
      log2foldchange > 1 ~ "red",      # Threshold for 'red'
      TRUE ~ "black"                   # Default to 'black'
    )
  )
  ```

## Saving Plots

- **Image Dimensions**: Modify the `ggsave()` function to adjust dimensions of saved plots by changing the `width` and `height` parameters.
  ```r
  ggsave(file.path(EXPERIMENTAL_FOLDER, "line_plot.png"), plot = line_plot, width = 8, height = 6)  # Adjust 'width' and 'height'
  ggsave(file.path(EXPERIMENTAL_FOLDER, "volcano_plot.png"), plot = volcano_plot, width = 8, height = 6)  # Adjust 'width' and 'height'
  ```

By following these instructions, you can effectively customize the plots in your R script to better fit your analytical needs and presentation requirements.

![image](https://github.com/user-attachments/assets/0179be74-56a4-4577-8207-5dbf1be378d2)

![volcano_plot](https://github.com/user-attachments/assets/d70d1959-62a4-426c-a840-0de4fa5e4c46)


