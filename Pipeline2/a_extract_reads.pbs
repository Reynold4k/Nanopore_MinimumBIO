#!/bin/bash

# Paths to directories and files
PARENT_DIR="/srv/scratch/z3546698/true/Small_Molecule/JQ1/T7MB-1/231104/R3/step3"
OUTPUT_DIR="/srv/scratch/z3546698/true/Small_Molecule/JQ1/T7MB-1/231104/potential_hit"
ANN_FILE="/srv/scratch/z3546698/true/reference/Homo_sapiens.GRCh38.110.gtf"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Extract all locations for hits from the GTF file and prefix with 'chr'
Hit_LOCATIONS=$(awk -v gene_name="BRD4" '$3 == "exon" && $0 ~ gene_name {print "chr"$1":"$4"-"$5}' "$ANN_FILE")

# Exit if no locations are found
if [ -z "$Hit_LOCATIONS" ]; then
  echo "No potential hits locations found in the annotation file."
  exit 1
fi

# Loop over each sorted BAM file
cd "$PARENT_DIR"
for bam_file in *_trimmed_sorted.bam; do
  # Check and create index if needed
  if [ ! -f "$bam_file.bai" ]; then
    echo "Indexing $bam_file..."
    samtools index "$bam_file"
  fi

  OUTPUT_BAM="${OUTPUT_DIR}/Hit_${bam_file}"
  echo "Extracting Potential hits from $bam_file to $OUTPUT_BAM..."

  # Extract reads for all potential hits locations
  samtools view -b "$bam_file" $Hit_LOCATIONS > "$OUTPUT_BAM"
done

echo "Extraction complete."
