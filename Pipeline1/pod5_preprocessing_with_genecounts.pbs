#!/bin/bash

# This script walks through the bioinformatics pipeline, starting from basecalling with Dorado, conducting quality control with NanoPlot,
# performing alignment, sorting, marking duplicates, and finally generating gene counts with featureCounts. Each step is crucial for
# ensuring high-quality, analyzable data from raw experimental results.

# What you need to modify:
# 1.Change the line 13 BASE_DIR to your/pod5/path
# 2.Change the line 14 REFERENCE to your/reference/path
# 3.Change the line 15 ANNOTATION to your/ANNOTATION/path

# Set the base directory and important paths
BASE_DIR="/srv/scratch/z3546698/true/pod5/"
REFERENCE="/srv/scratch/z3546698/true/reference/hg38.fa"
ANNOTATION="/srv/scratch/z3546698/true/reference/Homo_sapiens.GRCh38.110.gtf"
OUTPUT_DIR="$BASE_DIR/output"

# Load necessary modules on katana
module load dorado/0.7.3
module load porechop/0.2.4
module load nanoplot/1.43.0
module load bwa/0.7.17
module load samtools/1.20
module load subread/2.0.2
module load seqkit  # Ensure seqkit is available

# Create necessary directories
mkdir -p "$OUTPUT_DIR/fastq_pass"
mkdir -p "$OUTPUT_DIR/fastq_fail"
mkdir -p "$OUTPUT_DIR/step1"
mkdir -p "$OUTPUT_DIR/quality_control"
mkdir -p "$OUTPUT_DIR/step2"
mkdir -p "$OUTPUT_DIR/step3"
mkdir -p "$OUTPUT_DIR/step4"
mkdir -p "$OUTPUT_DIR/reports"

# Step 1: Merge .pod5 files and convert them to .fastq using Dorado
echo "Merging .pod5 files and converting them to .fastq..."
for STATE in pod5_pass pod5_fail; do

    # Directory to store the merged .pod5 files
    merged_pod5_dir="$OUTPUT_DIR/merged_$STATE"
    mkdir -p "$merged_pod5_dir"

    find "$BASE_DIR/$STATE" -type d -mindepth 1 -maxdepth 1 | while read -r barcode_dir; do
        merged_pod5_file="$merged_pod5_dir/$(basename "$barcode_dir").pod5"
        
        echo "Merging pod5 files in $barcode_dir..."
        # Assuming `npmerge` is used for merging .pod5 files
        npmerge "$barcode_dir"/*.pod5 -o "$merged_pod5_file"
        
        dorado_out="$OUTPUT_DIR/fastq_${STATE#pod5_}/$(basename "$barcode_dir")"
        mkdir -p "$dorado_out"
        
        # Basecalling with Dorado
        dorado basecaller --device "cuda" "$merged_pod5_file" "$dorado_out"
        
        echo "Converted: $merged_pod5_file -> $dorado_out"

        # Process with seqkit
        echo "Processing with seqkit for $dorado_out..."
        find "$dorado_out" -type f -name "*.fastq" | while read -r fastq_file; do
            seqprocessed_file="${fastq_file%.fastq}_processed.fastq"
            seqkit replace -p "^((?:.*?\n){3}).*?(GATCCGAATTC[ACGTN].*$)(\n.*)" -r '$1$2$3' -o "$seqprocessed_file" "$fastq_file"
            echo "Processed FASTQ file: $seqprocessed_file"
            mv "$seqprocessed_file" "$fastq_file"  # Overwrite the original with the processed file
        done
        
        # Optionally remove merged .pod5 files after conversion
        rm "$merged_pod5_file"
    done
done

# Step 2: Trimming fastq files using Porechop from fastq_pass directory
echo "Trimming fastq files using Porechop..."
for STATE in fastq_pass fastq_fail; do
    find "$OUTPUT_DIR/$STATE" -type f -name "*.fastq" | while read -r fastq_file; do
        step1_dir="$OUTPUT_DIR/step1/$(basename "$(dirname "$fastq_file")")"
        mkdir -p "$step1_dir"
        output_file="$step1_dir/$(basename "${fastq_file%.fastq}_trimmed.fastq.gz")"
        porechop -i "$fastq_file" -o "$output_file"
        echo "Trimmed: $fastq_file -> $output_file"
    done
done

# Step 3: Quality Control via NanoPlot
echo "Quality Control in progress..."
find "$OUTPUT_DIR/step1" -type f -name "*_trimmed.fastq.gz" | while read -r trimmed_file; do
    quality_control_dir="$OUTPUT_DIR/quality_control/$(basename "$(dirname "$trimmed_file")")"
    mkdir -p "$quality_control_dir"
    NanoPlot --fastq "$trimmed_file" --outdir "$quality_control_dir"
done

echo "Quality Control finished..."

# Step 4: Alignment and BAM file generation
echo "Alignment and BAM file generation in progress..."
find "$OUTPUT_DIR/step1" -type f -name "*_trimmed.fastq.gz" | while read -r trimmed_file; do
    dir="$OUTPUT_DIR/step2/$(basename "$(dirname "$trimmed_file")")"
    mkdir -p "$dir"
    basename=$(basename "$trimmed_file" .fastq.gz)
    output_bam="$dir/${basename}.bam"
    bwa mem "$REFERENCE" "$trimmed_file" | samtools view -Sb - > "$output_bam"
    echo "BAM file generated: $output_bam"
done

# Step 5: Sort and mark duplicates in BAM files
echo "Sorting and marking duplicates in BAM files..."
find "$OUTPUT_DIR/step2" -type f -name "*.bam" | while read -r bam_file; do
    dir="$OUTPUT_DIR/step3/$(basename "$(dirname "$bam_file")")"
    mkdir -p "$dir"
    sorted_bam="$dir/$(basename "${bam_file%.bam}_sorted.bam")"
    marked_bam="$dir/$(basename "${bam_file%.bam}_marked.bam")"

    # Sort BAM files
    samtools sort -o "$sorted_bam" "$bam_file"
    echo "Sorted BAM file: $sorted_bam"

    # Mark duplicates
    samtools markdup "$sorted_bam" "$marked_bam"
    echo "Marked duplicates in BAM file: $marked_bam"

    # Clean up intermediate files if needed
    rm "$sorted_bam"
done

# Step 6: Generate alignment statistics and reports
echo "Generating alignment statistics and reports..."
find "$OUTPUT_DIR/step3" -type f -name "*_marked.bam" | while read -r marked_bam; do
    dir="$OUTPUT_DIR/reports/$(basename "$(dirname "$marked_bam")")"
    mkdir -p "$dir"

    # Produce alignment statistics using samtools
    samtools flagstat "$marked_bam" > "$dir/$(basename "${marked_bam%.bam}_flagstat.txt")"
    echo "Alignment statistics report generated for: $marked_bam"
done

# Step 7: Feature counting with featureCounts
echo "Feature counting with featureCounts..."
find "$OUTPUT_DIR/step3" -type f -name "*_marked.bam" | while read -r marked_bam; do
    dir="$OUTPUT_DIR/step4/$(basename "$(dirname "$marked_bam")")"
    mkdir -p "$dir"
    
    output_counts="${dir}/$(basename "${marked_bam%.bam}_counts.txt")"
    featureCounts -a "$ANNOTATION" -o "$output_counts" -T 4 "$marked_bam"
    
    echo "Feature counts generated in: $output_counts"
done

echo "Pipeline completed. Please check the output directories for results and reports!"
