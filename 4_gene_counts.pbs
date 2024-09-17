#!/bin/bash

#Identify annotation file path
annotation="/mnt/d/hg38/hg38.ensGene.gtf"

FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders

# walk through all sorted and marked bam files
find "$FOLDER" -mindepth 1 -maxdepth 1 -type d | while read -r subfolder; do
    # extract folder names for the prefix of gene count matrix
    folder_name=$(basename "$subfolder")

    bam_files=$(find "$subfolder" -type f -name "*_sorted_marked.bam" -printf "%p ")

    # generate gene count matrix individually via featureCounts
    if [ -n "$bam_files" ]; then
        output_counts="$subfolder/${folder_name}_expression_counts.txt"

        # Perform feature counts
        featureCounts -a "$annotation" -o "$output_counts" -T 4 $bam_files

        echo "Feature counts for $folder_name are ready in: $output_counts"
    fi
done


