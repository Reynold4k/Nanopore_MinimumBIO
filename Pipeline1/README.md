

```bash
           " __  __ _____ _   _ _____ __  __ _    _ __  __ ____ _____ ____  "
           "|  \/  |_   _| \ | |_   _|  \/  | |  | |  \/  |  _ \_   _/ __ \ "
           "| \  / | | | |  \| | | | | \  / | |  | | \  / | |_) || || |  | |"
           "| |\/| | | | |     | | | | |\/| | |  | | |\/| |  _ < | || |  | |"
           "| |  | |_| |_| |\  |_| |_| |  | | |__| | |  | | |_) || || |__| |"
           "|_|  |_|_____|_| \_|_____|_|  |_|\____/|_|  |_|____/_____\____/ "
```

# Pipeline1

This document describes the Pipeline1 processing pipeline, which automates the analysis of sequencing data. It covers tasks like trimming, quality control, alignment, sorting, marking duplicates, and generating gene counts using various bioinformatics tools. The script is meant to be run on compatible systems, specifically on Linux environments or platforms like UNSW's Katana. This guide will help users customize the script as needed and understand where outputs are generated.

Important!! When running the processing and analysis scripts, you're supposed to change the working directorys first, and follow the step-to-step instructions.

## Features

1. **Trimming**: Utilizes Porechop for removing adapters from sequencing reads.
2. **Quality Control**: Leverages NanoPlot for assessing read quality.
3. **Alignment**: Aligns reads to a reference genome using BWA.
4. **Sorting and Deduplication**: Sorts and marks duplicate reads using Samtools.
5. **Gene Feature Counting**: Uses featureCounts for quantifying gene expression.

## Prerequisites

Ensure your system meets the following requirements before running the pipeline:

- **Operating System**: Linux or compatible environment.

## Rclone Example

To use a specific version of Rclone in your environment, you can list available versions and load the desired one as shown below:

## If you're using Katana UNSW, you can book an interactive CPU portal and then run the following command

```bash
qsub -I -l select=1:ncpus=16:mem=128gb -l walltime=04:00:00
```


### Notes

- **Module Availability**: Ensure all these modules are available on your current working path.

- Check them by command module av "the module"
  
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
```bash
sudo apt-get install porechop bwa samtools
```

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

## Index your reference fasta files, this step must be conducted on your current portal to continue the next step
```bash

bwa index hg38.fa

wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.gtf.gz 
gunzip Homo_sapiens.GRCh38.112.gtf.gz
```

# For katana only
4. Mount your OneDrive directories to your katana scratch, please refer to the official guidances:
https://docs.restech.unsw.edu.au/


### Usage Instruction

Set the input and reference data paths in the script:

Step1 Modify the preprocessing script(e.g. porechop_preprocessing.sh)
Do the same modification to the pod5 files processing scripts

For fastq files: porechop_preprocessing.pbs
For Pod5 files: pod5_preprocessing_with_genecounts.pbs


- **Modify the `FOLDER` variable**: Set it to the directory containing your sequencing files.

If you aren't sure about your location, you can print your location just by running:
```bash
pwd
```

```bash
FOLDER="the outout of your pwd working directory having the same structure as below"


For example,
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders
```

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

```bash 

bash ./porechop_preprocessing.sh

```

If run successfully, you'll see:
![image](https://github.com/user-attachments/assets/208572f3-bbff-4e30-acc2-107beacb8476)


### Clean Step: Edit the Cleaning Script and Execute

This script is designed to detect all expected output files generated from previous processing steps and, if all are confirmed to be present, it clears intermediate folders.
#### Execution:

1. Please edit the folder path with the same stratefy above in the step1,

For example,
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421" 

2. Save the script as `clear_intermediate.sh`.
3. Ensure the script has execution permissions by running:
   ```bash
   chmod +x clear_intermediate.sh

   ```
4. Execute the cleaning script via:
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

module load r/4.4.0

```

#Then edit R scripts in the gtf path:
```r
#Path to your newly generated Routput folder
plot_base_dir <- "/srv/scratch/z3546698/true/Routput"

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
# To find a GTF for Danio rerio (zebrafish), browse to:
# http://hgdownload.soe.ucsc.edu/goldenPath/danRer11/bigZips/genes/
# Example path for Danio rerio: "/path/to/your/directory/danRer11.refGene.gtf.gz"

# Make sure to decompress (.gz) files if needed and adjust the path accordingly.
```


Then run through the whole R scripts and check the result plots.
```bash

Rscript script.R
```

If running successfully, you'll see:

![image](https://github.com/user-attachments/assets/f19b438f-63be-456b-b67c-76ca02b3641a)









