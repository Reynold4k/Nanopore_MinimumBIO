#!/bin/bash

# Function to display a spinning progress indicator
spinner() {
    local pid=$!
    local delay=0.1
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep "$pid")" ]; do
        local temp=${spinstr#?}
        printf " [%c]  \b\b\b\b\b" "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
    done
    printf "    \b\b\b\b"
}

# Display the logo from logo.sh
source ./logo.sh

# Suppress function output using redirection
run_command() {
    ("$@") >/dev/null 2>&1 &
    spinner
}

# Step 1: Trimming and Quality Control using Porechop and NanoPlot
echo "Starting Step 1: Trimming and Quality Control..."
FOLDER="/mnt/d/Small_Molecule/Biotin/T7MB-2/240421"  # Path to input experiments folders

# Go to the working folder or print an error
cd "$FOLDER" || { echo "Folder not found"; exit 1; }

# Trimming fastq files using Porechop
echo "Trimming starts using Porechop......"

# Use find to locate all .fastq and .fastq.gz files in the directory and subdirectories
find . -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | while read -r file; do
    if [ -f "$file" ]; then
        # Generate output file name in the same directory
        output_file="${file%.fastq*}_trimmed.fastq.gz"
        # Run Porechop command
        run_command porechop -i "$file" -o "$output_file"
        echo "Trimmed: $file -> $output_file"
    fi
done

# Conducting Quality Control via NanoPlot
echo "Trimming finished......."
echo "Quality Control in progress......."

# Find and run NanoPlot on all trimmed files
find . -type f -name "*_trimmed.fastq.gz" | while read -r trimmed_file; do
    run_command NanoPlot --fastq "$trimmed_file" --outdir "${trimmed_file%.fastq.gz}_nanoplot"
done

echo "Step1 completed. Please check the sequencing quality report!"


# Step 2: Converting and Alignment
echo "Starting Step 2: Converting and Alignment..."
#!/bin/bash

# Identify folder path and reference path
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


echo "Step 2 completed."


# Step 3: Sorting and Marking Duplicates
echo "Starting Step 3: Sorting and Marking Duplicates..."
find "$FOLDER" -type f -name "*.bam" | while read -r bam_file; do
    dir=$(dirname "$bam_file")
    basename=$(basename "$bam_file" .bam)
    sorted_bam="$dir/${basename}_sorted.bam"
    marked_bam="$dir/${basename}_sorted_marked.bam"
    run_command samtools sort "$bam_file" -o "$sorted_bam"
    run_command samtools markdup "$sorted_bam" "$marked_bam"
    run_command samtools index "$marked_bam"
    rm "$sorted_bam"
    rm "$bam_file"
done
echo "Step 3 completed."

# Step 4: Generating Gene Counts and Annotation
echo "Starting Step 4: Gene Counts and Annotation..."
annotation="/mnt/d/hg38/hg38.ensGene.gtf"
find "$FOLDER" -mindepth 1 -maxdepth 1 -type d | while read -r subfolder; do
    folder_name=$(basename "$subfolder")
    bam_files=$(find "$subfolder" -type f -name "*_sorted_marked.bam" -printf "%p ")
    if [ -n "$bam_files" ]; then
        output_counts="$subfolder/${folder_name}_expression_counts.txt"
        run_command featureCounts -a "$annotation" -o "$output_counts" -T 4 $bam_files
    fi
done
echo "Step 4 completed."

# Notify user process is complete
echo "All steps completed"
source ./logo.sh
