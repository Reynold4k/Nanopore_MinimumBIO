#!/bin/bash

# Identify root path
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders

# find all bam files generated before and sort
find "$FOLDER" -type f -name "*.bam" | while read -r bam_file; do
    # acquire bam file path
    dir=$(dirname "$bam_file")
    
    # acquire filenames prefix
    basename=$(basename "$bam_file" .bam)

    
    sorted_bam="$dir/${basename}_sorted.bam"

    # sort bam files via samtools
    samtools sort "$bam_file" -o "$sorted_bam"

    echo "Sorted BAM file ready: $sorted_bam"

    # mark duplicates via samtools
    marked_bam="$dir/${basename}_sorted_marked.bam"
    samtools markdup "$sorted_bam" "$marked_bam"

    echo "Marked BAM file ready: $marked_bam"

    # create index files for all processed bam files
    samtools index "$marked_bam"

    echo "Indexed BAM file ready: $marked_bam"

    # rm all intermediate files to keep clear
    rm "$sorted_bam"
    rm "$bam_file"
done

