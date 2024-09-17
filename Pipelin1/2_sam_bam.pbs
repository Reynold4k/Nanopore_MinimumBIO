#!/bin/bash

# Identify folder path and reference path
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders
REFERENCE="/mnt/d/hg38/hg38.fa"

# Walk through each relevant FASTQ file for comparison and generate the BAM file
find "$FOLDER" -type f -name "*_trimmed.fastq.gz" | while read -r trimmed_file; do
    # trimmed file path
    dir=$(dirname "$trimmed_file")
    
    # acquiring filenames prefix
    basename=$(basename "$trimmed_file" .fastq.gz)

    # output bam files path
    output_bam="$dir/${basename}.bam"

    # conduct alignment and generate
    bwa mem "$REFERENCE" "$trimmed_file" | samtools view -Sb - > "$output_bam"

    echo "BAM file generated: $output_bam"
done

