#!/bin/bash

# Overview of the Bioinformatics Data Processing Pipeline

#This script is designed to process high-throughput sequencing data by executing a series of automated data processing steps on specified folders containing FASTQ files. The main functionalities include:
#1. **Merging and Trimming FASTQ Files**: It merges all FASTQ files within each sample directory into a single file and trims the merged file using Porechop to remove low-quality sequences and adapter sequences.
#2. **Quality Control**: Utilizes NanoPlot to assess the quality of the trimmed sequences, ensuring the reliability of the data.
#3. **Alignment and Coverage Analysis**: Leverages Minimap2 to align the trimmed sequences to a reference genome, generating BAM files and calculating read coverage, which is then output to `coverage.txt`.
#4. **Position Extraction and BED File Generation**: Extracts relevant gene positions from the coverage data and generates a BED file for easier visualization and downstream analysis.
#5. **Sequence Extraction**: Based on the generated BED file, extracts specified regions from the reference sequences and saves them in FASTA format for further analysis and applications.

#This script provides an automated workflow for bioinformatics research, significantly increasing data processing efficiency while ensuring accuracy and consistency in the results.


# What you what to modify:
# 1.Change the line 22 FASTQ_FOLDER to your/fastq/path
# 2.Change the line 27 REFERENCE to your/reference/path
# 3.Change the line 30 PROTEIN_DB to your/generated/protein/blast/database



FASTQ_FOLDER=(
    "/mnt/d/A_FKBP1B/WDB001/TON/241004"
)


REFERENCE="/mnt/d/hg38/Ton.fa"

# Path to the pre-built human protein BLAST database
PROTEIN_DB="/mnt/d/human_proteome_db"

# Define the pattern to match and replace
PATTERN="^((?:.*?\n){3}).*?(GATCCGAATTC[ACGTN].*$)(\n.*)"

for folder in "${FASTQ_FOLDER[@]}"; do
    logbook="${folder}/processing_logbook.txt"
    echo "Processing started for folder: $folder" > "$logbook"

    find "$folder" -mindepth 1 -maxdepth 1 -type d -name "R*" | while IFS= read -r sample_dir; do
        {
            echo "Processing directory: $sample_dir"

            # Define paths for step1 and step2
            step1_dir="${sample_dir}/step1"
            step2_dir="${sample_dir}/step2"

            # Create directories if they don't exist
            mkdir -p "$step1_dir"
            mkdir -p "$step2_dir"

            # Merge FASTQ files
            merged_file="${step1_dir}/all_sequences.fastq.gz"
            find "$sample_dir" -type f -name "*.fastq.gz" -print0 | xargs -0 cat > "$merged_file"

            if [ -s "$merged_file" ]; then
                echo "Files successfully merged into $merged_file."

                # Pre-process sequences with seqkit
                seqprocessed_file="${step1_dir}/all_sequences_processed.fastq.gz"
                seqkit replace -p "$PATTERN" -r '$1$2$3' -o "$seqprocessed_file" "$merged_file"
                echo "Sequences processed by seqkit and output to $seqprocessed_file."

                # Remove the original FASTQ files after processing, excluding files starting with 'all'
                find "$sample_dir" -type f -name "*.fastq.gz" ! -name "all*.fastq.gz" -delete
                echo "Original FASTQ files deleted from directory: $sample_dir (excluding 'all' prefixed files)"

                # Quality control using NanoPlot
                quality_control_dir="${sample_dir}/quality_control"
                mkdir -p "$quality_control_dir"
                NanoPlot --fastq "$seqprocessed_file" --outdir "${quality_control_dir}/"

                # Create a directory for BLASTp results
                blast_output_dir="${sample_dir}/blastp_results"
                mkdir -p "$blast_output_dir"

                # Translate FASTQ to protein sequences in a FASTA format
                protein_file="${blast_output_dir}/all_sequences_translated.fasta"
                zcat "$seqprocessed_file" | seqtk seq -A - | seqkit translate -o "$protein_file"
                echo "Translated sequences saved to $protein_file."

                # Run BLASTp
                output_blast="${blast_output_dir}/$(basename "${seqprocessed_file%.fastq.gz}_blastp.txt")"
                blastp -db "$PROTEIN_DB" -query "$protein_file" -out "$output_blast" -outfmt 6 -evalue 1e-3
                echo "BLASTp output saved to $output_blast."

                # Extract matching sequence IDs
                awk '{print $1}' "$output_blast" | sort | uniq > "${blast_output_dir}/matched_ids.txt"
                seqtk subseq "$seqprocessed_file" "${blast_output_dir}/matched_ids.txt" > "${step1_dir}/all_filtered_sequences.fastq.gz"
                echo "Filtered sequences saved to ${step1_dir}/all_filtered_sequences.fastq.gz if applicable."

                # Align reads to reference genome with Minimap2 and create BAM file for filtered sequences
                bam_file="${step2_dir}/aligned.bam"
                minimap2 -ax map-ont "$REFERENCE" "${step1_dir}/all_filtered_sequences.fastq.gz" | samtools view -Sb - > "$bam_file"

                # Sort BAM file
                sorted_bam_file="${step2_dir}/aligned_sorted.bam"
                samtools sort -o "$sorted_bam_file" "$bam_file"

                # Calculate coverage
                coverage_file="${step2_dir}/coverage.txt"
                samtools depth "$sorted_bam_file" > "$coverage_file"

                echo "Finished processing directory: $sample_dir. Results saved in step1 and step2."
            else
                echo "Merging failed for directory: $sample_dir. Skipping processing."
            fi

            # Handle extraction for each processed directory
            echo "Processing directory for extraction: $sample_dir"

            # Paths for necessary files
            coverage_file="${step2_dir}/coverage.txt"

            # Check if the coverage file exists before proceeding
            if [ -f "$coverage_file" ]; then
                echo "Preparing positions in new range format..."

                # Prepare ranges from coverage and generate positions.txt with coverage values
                awk '{
                    key = $1;  # Gene name
                    position = $2;  # Start position
                    coverage = $3;  # Coverage value

                    # Store the start and end positions in an array
                    if (key in positions) {
                        if (position < positions[key][1]) {
                            positions[key][1] = position;  # Update minimum start
                        }
                        positions[key][2] = position;  # Update to current position (end)
                        positions[key][3] += coverage;  # Accumulate coverage
                    } else {
                        positions[key][1] = position;  # Initialize min start
                        positions[key][2] = position;  # Initialize end
                        positions[key][3] = coverage;  # Initialize coverage
                    }
                } END {
                    for (gene in positions) {
                        # Print the gene with start-end range and total coverage
                        print gene "\t" positions[gene][1] "\t" positions[gene][2] "\t" positions[gene][3];
                    }
                }' "$coverage_file" > "${step2_dir}/positions.txt"

                # Sort to get the top 1000 genes by coverage
                sort -k4,4nr "${step2_dir}/positions.txt" | head -n 1000 > "${step2_dir}/top_1000_positions.txt"

                # Now create BED files from the top 1000 genes
                visualization_dir="${step2_dir}/visualization"
                mkdir -p "$visualization_dir"

                # Initialize a counter
                counter=0
                max_bed_files=1000

                while IFS=$'\t' read -r gene_id start end total_coverage; do
                    if [ "$counter" -ge "$max_bed_files" ]; then
                        echo "Reached maximum of $max_bed_files BED files. Stopping."
                        break
                    fi

                    # Generate the individual BED file name
                    bed_file="${visualization_dir}/${gene_id}.bed"

                    # Write the gene position to the new BED file
                    echo -e "$gene_id\t$start\t$end" > "$bed_file"
                    echo "Created BED file for $gene_id at $bed_file"

                    counter=$((counter + 1))
                done < "${step2_dir}/top_1000_positions.txt"
            else
                echo "Necessary coverage file not found in ${step2_dir}. Skipping directory."
            fi
        } >> "$logbook" 2>&1  # Redirect both stdout and stderr to the logbook
    done
done

echo "All directories processed." | tee -a "$logbook"
