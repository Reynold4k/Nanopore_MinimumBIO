# Pipeline 2

Once you have executed the initial analysis with your gene expression counts and identified potential hits, it’s time to further investigate them using AlphaFold for protein structure prediction. Below are detailed instructions for running the script while ensuring the relevant paths are properly set.

### Setting Up Paths
## Create a new conda environment
```bash
conda create -n alphafold_env python=3.9
conda activate alphafold_env
```

## Install canu, if it is not installed:
```bash
wget https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.bz2
tar -xvjf canu-2.2.Linux-amd64.tar.bz2
export PATH=$PATH:path/to/your/canu-2.2/bin
```
## Step A: Extracting potential hits into a new folder

#Edit file path in the script:

```bash
#Path to your bam files
PARENT_DIR="/path/to/your/experiments/R3/step2"
#Path to the new folder: potential hits
OUTPUT_DIR="/path/to/your/experiments/potential_hit"
#Annotation files
ANN_FILE="path/to/your/reference/Homo_sapiens.GRCh38.110.gtf"

Hit_LOCATIONS=$(awk -v gene_name="your_interested potential hits name" '$3 == "exon" && $0 ~ gene_name {print "chr"$1":"$4"-"$5}' "$ANN_FILE")

```

```bash

bash ./a_extract_reads.pbs

```

If running without error, you will see:

![image](https://github.com/user-attachments/assets/6114dda2-d07d-456c-9ac2-c96fe67e32f4)


## Step B: Error Correction, trimming, unitigging, Assembly

#Edit file path in the script:

OUTPUT_DIR="/path/to/your/experiments/potential_hit"

```bash

bash ./b_correction_canu.pbs

```

It won't take a long time, but you need to wait till the contigs files are created!!!

You must wait until seeing the contigs file:
![image](https://github.com/user-attachments/assets/9c6e2bd7-c371-4376-84dd-e84f345fe0f4)


Then execute the translation script, but still need to change the path first:

contig_fasta_file = "/path/to/your/corrected.contigs.fasta"
output_protein_file = "/path/to/your/protein_sequences_long.fasta"

```bash
python3 translation.py
```

If succeed, you will see:
![image](https://github.com/user-attachments/assets/e30b575e-6021-4724-9a16-80a87febc39c)

Make a blast online to check whether it is still the potential hits that you extracted in the Step A via: 

https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch


## Step C Finding overlapped region where could be indicated as the drug binding sites








