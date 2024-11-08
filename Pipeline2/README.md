## Pipeline 2
![image](https://github.com/user-attachments/assets/6d290d16-b830-485e-bf54-6cb27c42ef00)
![image](https://github.com/user-attachments/assets/57577119-5bcc-4c78-a291-996995d52a68)


Once you have executed the initial analysis with your gene expression counts and identified potential hits, it’s time to further investigate them using AlphaFold for protein structure prediction. Below are detailed instructions for running the script while ensuring the relevant paths are properly set.

### Setting Up 

## Create a new conda environment
```bash
conda create -n alphafold_env python=3.9

conda activate alphafold_env
```

This will download the latest version (v4) of the human proteome subset. Be sure you have enough space on your PC:
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
## Install required packages via:

```bash
#biopython
pip install biopython

#Canu for the error correction
apt install canu

#Bedtools2
apt-get install bedtools


```


## A quick implementation:
### (If you're using a personalized library, please skip PartA and PartB and use Pipeline2_1.5.sh and then Pipeline2_1_PartB_pymol.sh, and you may find your bed files in the visualization folder)

### For example:/srv/scratch/z3546698/tutorial/Bait_Glue/VHL/MB015/TON/231124/R1/step2/visualization

# PartA

# Pipeline 1
## After setting up the environment and libraries in the stepD, you need to change several directories setting in the script: pipeline2_1_PartA.sh

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

bash ./pipeline2_1_PartA.sh
```

## Then some operations are also required for part B:

For the pipeline2_1_PartB.sh, you need to change:

```bash

GENOME_FASTA="Path to reference/hg38.fa"  # Reference FASTA file

#choose one gene file each time only to avoid errors:
BED_FILE="Path to your selected gene bed file/FKBP1C_Hit_all_trimmed_sorted_merged.bed"  # BED file with high coverage regions

bash ./pipeline2_1_PartB.sh
```
# Pipeline 1.5

## What you need to modify:
1.Change the line 15 GENOME_FASTA to your/reference/path
```bash
GENOME_FASTA="/mnt/c/Users/70921/OneDrive/Desktop/reference/T7-Pep_Ref_93nt.fasta" # Reference FASTA file
```
2.Change the line 16 BED_FILE to your/bed/path
```bash
# Define files and directories
BED_FILE="/mnt/d/Bait_Glue/VHL/MB012/TON/230827/R1/visualization/NP_006156.2_2.bed"  # BED file
```

# PartB
# Then both pipeline1 and pipeline1.5 are supposed to take some operations are also required for part B_pymol:

For the pipeline2_1_PartB_pymol.sh, you need to change:

```bash

VISUAL_DIR="path to your potential_hit/visualization"  # Directory for visualizations

#Or if you're using your personalized library, just specify it to your experimental folder to view the result faster

#For example:
VISUAL_DIR="/srv/scratch/z3546698/tutorial/Bait_Glue/VHL/MB015/TON/231124"  # Directory for visualizations



GENOME_FASTA="Path to reference/hg38.fa"  # Reference FASTA file

# Search the gene at uniprot website to check the right uniprot Id for each gene and edit:
UNIPROT_ID="search at unitprot website"
https://www.uniprot.org/

#Path to your prepared alphafold database library in Step D!!:

PDB_PATH="/srv/scratch/z3546698/true/alphafold/database/UP000005640_9606_HUMAN_v4/AF-${UNIPROT_ID}-F1-model_v4.pdb"

bash ./pipeline2_1_PartB_pymol.sh

```




