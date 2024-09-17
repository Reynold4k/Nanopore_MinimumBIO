

```bash
           " __  __ _____ _   _ _____ __  __ _    _ __  __ ____ _____ ____  "
           "|  \/  |_   _| \ | |_   _|  \/  | |  | |  \/  |  _ \_   _/ __ \ "
           "| \  / | | | |  \| | | | | \  / | |  | | \  / | |_) || || |  | |"
           "| |\/| | | | |     | | | | |\/| | |  | | |\/| |  _ < | || |  | |"
           "| |  | |_| |_| |\  |_| |_| |  | | |__| | |  | | |_) || || |__| |"
           "|_|  |_|_____|_| \_|_____|_|  |_|\____/|_|  |_|____/_____\____/ "
```


# Comprehensive Pipeline for Trimming, Quality Control, Alignment, and Analysis

**Author**: Chen Zhu  
**Email**: [z3546698@ad.unsw.edu.au]
**Date**: 04/Sep/2024

## Pipeline1

This pipeline automates the process of trimming, quality control, alignment, sorting, marking duplicates (if applicable), and generating gene counts from sequencing data. It uses several bioinformatics tools to achieve efficient processing and analysis.

## Features

1. **Trimming**: Utilizes Porechop for removing adapters from sequencing reads.
2. **Quality Control**: Leverages NanoPlot for assessing read quality.
3. **Alignment**: Aligns reads to a reference genome using BWA.
4. **Sorting and Deduplication**: Sorts and marks duplicate reads using Samtools.
5. **Gene Feature Counting**: Uses featureCounts for quantifying gene expression.

## Prerequisites

Ensure your system meets the following requirements before running the pipeline:

- **Operating System**: Linux or compatible environment.

**Software requirements**:
## Software List

1. **Porechop**   (1_porechoptrim_fastqc.sh)
   - Tool for trimming Oxford Nanopore reads.
   - Module command: `module load porechop`

2. **Cutadapt**
   - Tool for trimming Oxford Nanopore reads.
   - Module command: `module load cutadapt`

3. **Dorado**     (dorado_preprocessing.sh)
   - Tool for trimming Oxford Nanopore reads.
   - Module command: `module load dorado`

4. **NanoPlot**   (1_porechoptrim_fastqc.sh)
   - Tool for quality control of nanopore reads.
   - Module command: `module load nanoplot`

5. **BWA**        (2_sam_bam.sh)
   - Aligns sequence reads to a reference genome.
   - Module command: `module load bwa`

6. **Samtools**   (2_sam_bam.sh, 3_sort_markdup.sh)
   - Utilities for manipulating alignments in the SAM format, including sorting and indexing.
   - Module command: `module load samtools`

7. **featureCounts**    (4_gene_counts.sh)
   - Part of the Subread package for counting reads to genomic features.
   - Module command: `module load subread`

For the script <porechop_preprocessing.sh>, please load 1,4,5,6,7
For the script <dorado_preprocessing.sh>, please load 3,4,5,6,7


## Rclone Example

To use a specific version of Rclone in your environment, you can list available versions and load the desired one as shown below:

## If you're using Katana UNSW, you can book an interactive CPU portal and then run the following command

qsub -I -l select=1:ncpus=16:mem=128gb -l walltime=04:00:00


```bash
# Display all Rclone versions available
module avail your_module

# Load the desired version(porechop, nanoplot, bwa, samtools, subread)
module load your_module/version_number

# It may also be the case that some modules are not included in katana system packages, then please follow the procedure below:

module load python/3.8.15  # please specify a python environment, for example python version 3.8.15

# Create and Activate virtual environment for running pip command (Highly Recommended)

python -m venv nanoplot_env
source nanoplot_env/bin/activate

# Please pip install the packages that you cannot find on katana
pip install NanoPlot, dorado, bwa, ...

# Using --help to check the availability of installed packages
NanoPlot --help
```


### Notes

- **Module Availability**: Ensure all these modules are available on your current working path.
  
- **Version Compatibility**: Always check compatibility of module versions with your data and analysis workflow.


**Reference Files**:
- Reference genome in FASTA format (e.g., `hg38.fa`).
- Gene annotation file in GTF format (e.g., `hg38.ensGene.gtf`).

## Installation Dependencies (If you're working on your local environment such as Ubuntu, if you're working with Katana, skip to usage instructions)



To install necessary bioinformatics tools, you might require package managers like `conda` or `apt-get`:

```bash
# Create working environment with conda
conda create -n bioenv porechop nanopolish bwa samtools subread
conda activate bioenv

```

# Installation with apt-get
sudo apt-get install porechop bwa samtools

## Setup Instructions

1. Clone the repository or copy the script into your working directory.
2. Ensure that the script and any additional scripts (`logo.sh`) are in the same directory.
3. Download hg38 reference and annotation files via:

```bash

mkdir -p /path/to/reference
cd /path/to/reference
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

```

# This step usually takes around an hour even on Katana
bwa index hg38.fa

wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz 
gunzip Homo_sapiens.GRCh38.112.gtf.gz


# For katana only
4. Mount your OneDrive directories to your katana scratch, please refer to the official guidances:
https://docs.restech.unsw.edu.au/


### Usage Instruction

Set the input and reference data paths in the script:

Step1 Modify the preprocessing script(e.g. porechop_preprocessing.sh)

For fastq files: porechop_preprocessing.sh
For Pod5 files: pod5_preprocessing.sh (Without gene counts)
For Pod5 files: pod5_preprocessing_with_genecounts.sh (With gene counts)

- **Modify the `FOLDER` variable**: Set it to the directory containing your sequencing files.

If you aren't sure about your location, you can print your location just by running:

pwd


FOLDER="the outout of your pwd working directory having the same structure as below"

For example,
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders


Step2 Update gene reference and annotation file path in the pre processing script

- **Update the `REFERENCE` and `annotation` paths**: Point these to your reference genome and annotation files, respectively.

For example,
REFERENCE="/mnt/d/hg38/hg38.fa"
annotation="/mnt/d/hg38/hg38.ensGene.gtf"

If you are using UNSW KATANA, the working directory format may differ

Ensure your directory structure allows the script to find input files and output directories correctly as described in the script comments.



Step3 Save the edited script and execute script

- **Execute the scripts for their corresponding data

For example, for nanopore fastq data, run command:

```bash 

./porechop_preprocessing.sh

```

## Trouble Shooting

If encounter errors or no responses when running the scripts, please run split scripts step by step.

The porechop_preprocessing.sh is consisted of 4 split scripts:

1_porechoptrim_fastqc.sh
2_sam_bam.sh
3_sort_markdup.sh
4_gene_counts.sh
clear_intermediate.sh


### Clean Step: Edit the Cleaning Script and Execute

This script is designed to detect all expected output files generated from previous processing steps and, if all are confirmed to be present, it clears intermediate folders (`step1` to `step4`) within the specified directory.

Please edit the folder path with the same stratefy above in the step1,

For example,
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421" 

After edited the folder path in the script, excute via:
```bash

./clear_intermediate.sh

```

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

#### Execution:
1. Save the script as `clear_intermediate.sh`.
2. Ensure the script has execution permissions by running:
   ```bash
   chmod +x clear_intermediate.sh

   ```


## Example File Path Setup

The script expects input files to be in the specified `FOLDER` directory and also looks for `.fastq` or `.fastq.gz` files for trimming and alignment.

Input folder path example and its supposed structure:
/mnt/d/Small_Molecule/Biotin/T7MB-2/240421
  ├── R0
  ├── R1
  └── ...
       ├── sample1.fastq.gz
       ├── sample2.fastq.gz
       └── ...

Output folder structure should be similar as below:

/mnt/d/Small_Molecule/Biotin/T7MB-2/240421
  ├── R0
  ├── R1
  └── ...
       ├── sample.fastq.gz
       ├── sample_trimmed.fastq.gz
       ├── sample_trimmed.bam
       ├── **FAY71653_pass_barcode20_8be234fc_b19da987_8_trimmed_nanoplot**
       ├── sample_trimmed_sorted_marked.bam
       ├── sample_trimmed_sorted_marked.bam.bai
       ├── sample_expression_counts.txt
       ├── sample_expression_counts.txt.summary
       └── ...

Output for Pod5 data:

Path/to/your/pod5/output
├── fastq_fail
    └── barcode....
├── fastq_pass
│   └── barcode....
├── quality_control
├── step1
├── step2
├── step3
├── step4
├── step5
└── step6

## Output

The pipeline produces the following outputs:

- **Trimmed Fastq Files**: Located in the same directory as input with `_trimmed.fastq.gz` suffix.
- **Quality Reports**: Generated by NanoPlot in directories named as `<sample>_nanoplot`.
- **BAM Files**: Alignments stored as `<sample>.bam` files.
- **Sorted and Deduplicated BAM Files**: Named `<sample>_sorted_marked.bam`.
- **Expression Counts**: Gene expression data saved as `<folder_name>_expression_counts.txt`.



## Analysis of generated gene counts matrix: Updating 3 Paths in R and Bash Scripts

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

#Then edit R scripts in the gtf path:
```r

gtf_file <- "/path/to/your/directory/Homo_sapiens.GRCh38.112.gtf.gz"
```


Then run through the whole R scripts and check the result plots.
```bash

Rscript script.R
```









