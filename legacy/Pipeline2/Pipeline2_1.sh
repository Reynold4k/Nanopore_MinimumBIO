#!/bin/bash

# This script automates the analysis and visualization of potential genomic hits from RNA sequencing data focused on specific genes,
# and subsequently extracts DNA sequences from specified genomic coordinates and translates them into protein sequences.

# Part 1: RNA-seq Analysis

# Paths to directories and files
PARENT_DIR="/mnt/d/A_FKBP1B/WDB001/YB/241004/R3/step2"
GENOME_FASTA="/mnt/d/hg38/hg38.fa"  # Path to the reference genome FASTA file
ANN_FILE="/mnt/d/hg38/Homo_sapiens.GRCh38.112.gtf"  # Path to the GTF file for gene annotations

OUTPUT_DIR="${PARENT_DIR%/*/*}/potential_hit"
VISUAL_DIR="${OUTPUT_DIR}/visualization"
FASTQ_DIR="${OUTPUT_DIR}/FASTQ"

# Ensure necessary directories exist
mkdir -p "$VISUAL_DIR" "$FASTQ_DIR" "$OUTPUT_DIR"

# Specify a list of target genes to analyze
GENES=("PHF6")

# Change to the PARENT_DIR directory
cd "$PARENT_DIR" || { echo "Failed to access directory: $PARENT_DIR"; exit 1; }

# Define the BAM file to process
bam_file="all_filtered_sequences_sorted.bam"

# Check if the BAM file exists
if [ ! -f "$bam_file" ]; then
  echo "Error: $bam_file does not exist in $PARENT_DIR."
  exit 1
fi

# Ensure the BAM file has an index; create it if missing
if [ ! -f "$bam_file.bai" ]; then
  echo "Indexing $bam_file..."
  samtools index "$bam_file"
fi

# Extract data based on locations found in the GTF file
for gene_name in "${GENES[@]}"; do
  # Adjust the chromosome naming if needed by adding 'chr' prefix
  Hit_LOCATIONS=$(awk -v gene="$gene_name" '$3 == "exon" && $0 ~ gene {print "chr"$1":"$4"-"$5}' "$ANN_FILE")

  if [ -z "$Hit_LOCATIONS" ]; then
    echo "No potential hit locations found for $gene_name in the annotation file."
    continue
  fi

  OUTPUT_BAM="${OUTPUT_DIR}/${gene_name}_Hit_${bam_file}"
  echo "Extracting potential hits for $gene_name from $bam_file to $OUTPUT_BAM..."
  samtools view -b "$bam_file" $Hit_LOCATIONS > "$OUTPUT_BAM"

  # Convert to BED format
  bed_file="${OUTPUT_BAM%.bam}.bed"
  bedtools bamtobed -i "$OUTPUT_BAM" > "$bed_file"

  # Sort the BED file
  sorted_bed_file="${bed_file%.bed}_sorted.bed"
  sort -k1,1 -k2,2n "$bed_file" > "$sorted_bed_file"

  echo "Bed file for $gene_name generated at $sorted_bed_file"
done

# Ensure the script uses the correct BED file for DNA extraction
BED_FILE="$sorted_bed_file"  # Use the latest sorted BED file

# Set output path for DNA sequences
dna_sequences_output="${OUTPUT_DIR}/high_coverage_sequences.fa"
translated_proteins_output="${OUTPUT_DIR}/translated_proteins.fa"

# Check the integrity of the BED_FILE
if [ ! -f "$BED_FILE" ]; then
  echo "Error: Expected BED file $BED_FILE not found."
  exit 1
fi

# Proceed to BED-based DNA extraction
bedtools getfasta -fi "$GENOME_FASTA" -bed "$BED_FILE" -fo "$dna_sequences_output" ||
{ echo "Error: Could not extract DNA sequences."; exit 1; }

echo "DNA sequences extracted to $dna_sequences_output"

# Translate DNA to protein using Python with Biopython
if [ -f "$dna_sequences_output" ]; then
  python3 - <<EOF
from Bio.Seq import Seq
from Bio import SeqIO

input_file = "$dna_sequences_output"
output_file = "$translated_proteins_output"

with open(output_file, "w") as f_out:
    for record in SeqIO.parse(input_file, "fasta"):
        dna_seq = Seq(str(record.seq))
        print(f"Processing sequence {record.id} of length {len(dna_seq)}")

        if len(dna_seq) % 3 != 0:
            pad_length = (3 - len(dna_seq) % 3) % 3
            dna_seq = dna_seq + 'N' * pad_length
            print(f"Padding sequence {record.id} to a length of {len(dna_seq)}")

        protein_seq = dna_seq.translate(to_stop=False)

        # Write the translated protein sequence to output
        f_out.write(f">{record.id}\n{str(protein_seq)}\n")

print(f"Protein sequences saved to {output_file}")
EOF

echo "Protein translation complete. Sequences saved to $translated_proteins_output."
else
  echo "Error: No DNA sequences found in $dna_sequences_output"
fi
