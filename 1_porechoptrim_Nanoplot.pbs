#!/bin/bash

# Identify working path and primers
FOLDER="/srv/scratch/z3546698/true/Bait_Glue/DCAF15/No Glue/T7MB-2/240504"

# Go to the working folder or print an error
cd "$FOLDER" || { echo "Folder not found"; exit 1; }

# Trimming fastq files using Porechop
echo "Trimming starts using Porechop......"

# Use find to locate all .fastq and .fastq.gz files in the directory and subdirectories
find . -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while read -r file; do
    if [ -f "$file" ]; then
        # Generate output file name in the same directory
        output_file="${file%.fastq*}_trimmed.fastq.gz"  # Change to retain the original folder structure
        # Make sure output name doesn’t conflict
        porechop -i "$file" -o "$output_file"  # Use Porechop command
        echo "Trimmed: $file -> $output_file"
    fi
done

# Conducting Quality Control via NanoPlot software
echo "Trimming finished......."
echo "Quality Control in progress......."

# Find and run FastQC on all trimmed files
find . -type f -name "*_trimmed.fastq.gz" | while read -r trimmed_file; do
    NanoPlot --fastq "$trimmed_file" --outdir "${trimmed_file%.fastq.gz}_nanoplot"
done

echo "All pre-processing done. Please check the sequencing quality report!"
