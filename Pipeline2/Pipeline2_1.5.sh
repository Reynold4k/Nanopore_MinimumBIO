#!/bin/bash

# This script automates the analysis and visualization of potential genomic hits from RNA sequencing data focused on specific genes.
# It is structured to process BAM files to extract regions of interest, convert these into FASTQ files, and perform a genome assembly
# using Canu for possible error correction. The process includes extracting hit regions from annotation files, converting data formats,
# and visualizing coverage. Overall, this pipeline is designed to facilitate gene-specific analyses within high-throughput datasets,
# allowing researchers to focus on particular genes of interest, manage raw data complexity, and enhance sequence accuracy with assembly strategies.

# What you need to modify:
# 1.Change the line 15 GENOME_FASTA to your/reference/path

# 2.Change the line 16 BED_FILE to your/bed/path

# Define files and directories
GENOME_FASTA="/mnt/c/Users/70921/OneDrive/桌面/reference/T7-Pep_Ref_93nt.fasta" # Reference FASTA file
BED_FILE="/mnt/d/Bait_Glue/VHL/MB012/TON/230827/R1/visualization/NP_006156.2_2.bed"  # BED file

# Extract directory path for the outputs
BED_DIR=$(dirname "$BED_FILE")

# Define output paths
dna_sequences_output="$BED_DIR/high_coverage_sequences.fa"
translated_proteins_output="$BED_DIR/translated_proteins.fa"

# Initialize output file
echo "" > "$dna_sequences_output"

# Read the BED file and extract sequences from the FASTA
while read -r line; do
    # Skip comments and empty lines
    [[ "$line" =~ ^# || -z "$line" ]] && continue

    # Read BED fields, assuming they are delimited by whitespace (could also be tab)
    chrom=$(echo $line | awk '{print $1}')    # 1: Refseq (ID)
    start=$(echo $line | awk '{print $2}')    # 2: start position
    end=$(echo $line | awk '{print $3}')      # 3: end position

    # Extract the relevant sequence from the FASTA file using a custom Python script
    python3 - <<EOF
from Bio import SeqIO

genome_fasta = "$GENOME_FASTA"
output_file = "$dna_sequences_output"

# Find the sequence we want to extract
chrom = "$chrom"
start = int($start)  # Convert to integer (0-based)
end = int($end)

# Read the genome FASTA and extract the requested sequence
with open(genome_fasta, "r") as fasta_file, open(output_file, "a") as f_out:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == chrom:
            # Extract the desired sequence (remember Python uses 0-based index)
            extracted_seq = record.seq[start:end]
            f_out.write(f">{record.id}:{start}-{end}\n{extracted_seq}\n")

EOF

done < "$BED_FILE"

echo "DNA sequences extracted to $dna_sequences_output"

python3 - <<EOF
from Bio.Seq import Seq
from Bio import SeqIO

input_file = "$dna_sequences_output"
output_file = "$translated_proteins_output"

with open(output_file, "w") as f_out:
    for record in SeqIO.parse(input_file, "fasta"):
        dna_seq = Seq(str(record.seq))

        # Print sequence length for debugging
        print(f"Processing sequence {record.id} of length {len(dna_seq)}")

        # Handling sequences not a multiple of three
        if len(dna_seq) % 3 != 0:
            # Pad with N's to make length a multiple of three
            pad_length = (3 - len(dna_seq) % 3) % 3
            dna_seq = dna_seq + 'N' * pad_length
            print(f"Padding sequence {record.id} to a length of {len(dna_seq)}")

        protein_seq = dna_seq.translate(to_stop=False)

        # Write the translated protein sequence to output
        f_out.write(f">{record.id}\n{str(protein_seq)}\n")

print(f"Protein sequences saved to {output_file}")
EOF
