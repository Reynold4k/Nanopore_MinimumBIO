## Pipeline 2

Once you have executed the initial analysis with your gene expression counts and identified potential hits, it’s time to further investigate them using AlphaFold for protein structure prediction. Below are detailed instructions for running the script while ensuring the relevant paths are properly set.

### Setting Up Paths
## Create a new conda environment
```bash
conda create -n alphafold_env python=3.9
conda activate alphafold_env
```

## A quick implementation:

## After setting up the environment and libraries in the stepD, you need to change several directories setting in the script: pipeline2-A.pbs

```bash
PARENT_DIR="/path/to/your/bam_files_directory"
# For example:
PARENT_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/R3/step2"

OUTPUT_DIR="/path/to/your/output_directory"
# For example:
OUTPUT_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/potential_hit"

ANN_FILE="/path/to/your/reference/Homo_sapiens.GRCh38.110.gtf"
# For example:

ANN_FILE="path to /Homo_sapiens.GRCh38.110.gtf"

#Gene names!! Check the spelling with NCBI to avoid gene name spellin errors:

# For example:
GENES=("FKBP1A" "FKBP1C") # Add your desired gene names here

# Create a new visualization path in your working directory:
VISUAL_DIR="path to your potential_hit/visualization"  # Directory for visualizations

bash ./pipeline2-A.pbs
```

Then some operations are also required for part B:

```bash
GENOME_FASTA="Path to reference/hg38.fa"  # Reference FASTA file

#choose one gene file each time only to avoid errors:
BED_FILE="Path to your selected gene bed file/FKBP1C_Hit_all_trimmed_sorted_merged.bed"  # BED file with high coverage regions

VISUAL_DIR="path to your potential_hit/visualization"  # Directory for visualizations

# Search the gene at uniprot website to check the right uniprot Id for each gene and edit:
UNIPROT_ID="Q5VVH2"

#Path to your prepared alphafold database library in Step D!!:

PDB_PATH="/srv/scratch/z3546698/true/alphafold/database/UP000005640_9606_HUMAN_v4/AF-${UNIPROT_ID}-F1-model_v4.pdb"

```


## Below are the instructions for manually running the scripts step by step, you don't need step A-C scripts for now, if you're keen on how it works, please contact me at my email.


## Step A: Extracting potential hits into a new folder

#Edit file path in the script:

# Extract Hits Script

This guide provides instructions to set up the environment and run a script to extract potential hits from BAM files for specified genes, using a reference GTF file. Before running the script, ensure you have the necessary software installed and the environment set up correctly.

## Prerequisites

- **Bash**: Ensure your system can execute bash scripts.
- **Samtools**: The script requires Samtools version 1.20. You can install it via package managers or build it from source.

### Installation Instructions for Samtools

1. **Using a Package Manager**:
   - On Ubuntu: `sudo apt-get install samtools`
   - On CentOS/RHEL: `sudo yum install samtools`

2. **Build from Source**:
   - Download the source code from [Samtools Releases](https://github.com/samtools/samtools/releases).
   - Follow the instructions in the `INSTALL` file provided with the source code.

3. **Loading via Module (HPC Systems)**:
   - If using an HPC system, load Samtools using the module system:
     ```bash
     module load samtools/1.20
     ```

## Script Configuration

Before running the script, update the following paths to match your file system's structure:

- **PARENT_DIR**: Path to the directory containing your sorted BAM files.

```bash
PARENT_DIR="/path/to/your/bam_files_directory"
# For example:
PARENT_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/R3/step2"

OUTPUT_DIR="/path/to/your/output_directory"
# For example:
OUTPUT_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/potential_hit"

ANN_FILE="/path/to/your/reference/Homo_sapiens.GRCh38.110.gtf"
# For example:

ANN_FILE="/srv/scratch/z3546698/true/reference/Homo_sapiens.GRCh38.110.gtf"

#Gene names!! Check the spelling with NCBI to avoid gene name spellin errors:

# For example:
GENES=("FKBP1A" "FKBP1C") # Add your desired gene names here


```

## Step B: Error Correction, trimming, unitigging, Assembly

#Edit file path in the script:

OUTPUT_DIR="/path/to/your/experiments/potential_hit"

```bash

bash ./b_correction_canu.pbs

```

It won't take a long time, but you need to wait till the contigs files are created!!!

You must wait until seeing the contigs file:
![image](https://github.com/user-attachments/assets/ffd066f3-b994-4549-9a62-f2993a58fa5c)

## Step C Finding overlapped region where could be indicated as the drug binding sites

To continue with step C, you still need to modify the directories and check the results before the next step,
This will ensure you enough confidence and no errors encountered:

```bash
# Define directories
OUTPUT_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/potential_hit"
VISUAL_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/potential_hit/visualization"

```



Sample Illustration:

![image](https://github.com/user-attachments/assets/a686583f-6922-49f9-978f-19294cd709c4)


## Step D Mapping your extracted gene reads segments to alphafold library


This will download the latest version (v4) of the human proteome subset. Be sure you have enough space (e.g. in the `scratch` directory in Katana):
```bash
# Navigate to the target directory
cd /srv/scratch/z3546698/true/alphafold/database/

# Download the AlphaFold database tar file if not already downloaded
wget -c https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar

# Extract the .tar file
tar -xvf UP000005640_9606_HUMAN_v4.tar

mkdir UP000005640_9606_HUMAN_v4

# Navigate into the extracted directory
cd UP000005640_9606_HUMAN_v4

# Remove .cif.gz files (if any exist)
rm *.cif.gz

# Unzip .pdb.gz files
gunzip *.gz
```

## Install Pymol

```bash

#Install Pymol if it is not available to load at your terminal:
conda create -n pymol_env python=3.8
conda activate pymol_env
conda install -c conda-forge -c schrodinger pymol-bundle

```

## Install JQ
```bash

wget https://github.com/jqlang/jq/releases/download/jq-1.7.1/jq-1.7.1.tar.gz
tar -xzvf jq-1.7.1.tar.gz
cd jq-1.7.1
autoreconf -i
./configure
make

#Test if it is available now, it will work even some parts go wrong during installation, just check it 
jq --version

```



