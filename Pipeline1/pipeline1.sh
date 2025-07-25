#!/bin/bash

# This script outlines a comprehensive bioinformatics pipeline designed to preprocess, analyze, and interpret high-throughput sequencing data. 
# The process begins with merging and trimming FASTQ files to remove adapters and clean up the sequence data, ensuring quality inputs for alignment. 
# Trimming is performed using Porechop, followed by Quality Control assessments with NanoPlot to verify data integrity and readiness for downstream analysis. 
# The next step involves aligning the cleaned reads to a reference genome using BWA, producing BAM files that map the reads appropriately within genomic coordinates.
# These BAM files are then sorted, and duplicates are marked to reduce biases incidental to PCR amplification, using tools such as Samtools. 
# Finally, featureCounts is used to quantify gene expression levels by measuring how many reads align to each genomic feature, such as a gene or exon, 
# providing key data for differential expression studies and other genomic investigations. 
# This pipeline integrates essential preprocessing steps to ensure that the data is of high quality and ready for comprehensive genomic analyses.

# What you need to modify:
# 1.Change the line 19 FOLDER to your/fastq/path
# 2.Change the line 25 REFERENCE to your/reference/path
# 3.Change the line 27 ANNOTATION to your/ANNOTATION/path
# 4.Change the line 32 PROTEIN_DB to your/makeblastdb/protein/database/path


FOLDERS=(
    "/mnt/d/A_FKBP1B/No_glue/YB/241004"
)

REFERENCE="/mnt/d/hg38/hg38.fa"
ANNOTATION="/mnt/d/hg38/Homo_sapiens.GRCh38.112.gtf"
PROTEIN_DB="/mnt/d/human_proteome_db"

for FOLDER in "${FOLDERS[@]}"; do
  echo "Merging and processing files in folder $FOLDER......"
  LOGFILE="$FOLDER/script_output.log"
  exec > >(tee -a "$LOGFILE") 2>&1  # Redirect stdout and stderr to the log file

  find "$FOLDER" -type d -name "R*" | while IFS= read -r dir; do
      echo "Processing directory: $dir"
      merged_file="$dir/all_sequences.fastq.gz"

      if [ ! -e "$merged_file" ]; then
          echo "Merging files..."
          temp_merged="${merged_file%.gz}"  # Create temp file without .gz extension
          > "$temp_merged"  # Clear or initialize the temp file

          find "$dir" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -print0 | while IFS= read -r -d '' file; do
              echo "Merging file: $file"
              if [[ "$file" == *.gz ]]; then
                  zcat "$file" >> "$temp_merged"
              else
                  cat "$file" >> "$temp_merged"
              fi
          done
          echo "Compressing merged file..."
          gzip "$temp_merged"  # Gzip the merged file
          echo "Files merged into $merged_file."
      else
          echo "Merged file $merged_file already exists. Skipping merging."
      fi

      seqprocessed_file="$dir/all_sequences_processed.fastq.gz"
      if [ ! -e "$seqprocessed_file" ]; then
          if cutadapt -j 16 -g "GATCCGAATTCN" -e 0.01 --times=3 --discard-untrimmed -o "$seqprocessed_file" "$merged_file"; then
              echo "Sequences processed and saved to $seqprocessed_file"
          else
              echo "Error: cutadapt failed to process sequences."
              continue  # Go to the next directory
          fi
      else
          echo "Processed file $seqprocessed_file already exists. Skipping processing."
      fi

      step1_dir="$dir/step1"
      mkdir -p "$step1_dir"
      trimmed_file="$step1_dir/all_trimmed.fastq.gz"

      if [ ! -e "$trimmed_file" ]; then
          if cp "$seqprocessed_file" "$trimmed_file"; then
              rm "$seqprocessed_file"
              echo "Trimmed file created: $trimmed_file"
          else
              echo "Error: Unable to copy processed file $seqprocessed_file."
          fi
      else
          echo "Trimmed file $trimmed_file already exists. Skipping copy."
      fi
  done

  echo "Quality Control in progress......."
  find "$FOLDER" -type f -path "*/step1/*_trimmed.fastq.gz" | while read -r trimmed_file; do
      quality_control_dir="$(dirname "$trimmed_file")/../quality_control"
      mkdir -p "$quality_control_dir"
      nanop_output_dir="$quality_control_dir/$(basename "${trimmed_file%.fastq.gz}_nanop")"

      if [ ! -d "$nanop_output_dir" ]; then
          NanoPlot --fastq "$trimmed_file" --outdir "$nanop_output_dir"
          echo "Quality control data generated in $nanop_output_dir"
      else
          echo "Quality control data already exists at $nanop_output_dir. Skipping quality control."
      fi
  done

  echo "BLASTp in-frame check in progress......."
  find "$FOLDER" -type f -path "*/step1/*_trimmed.fastq.gz" | while read -r trimmed_file; do
      dir=$(dirname "$trimmed_file")
      blast_output_dir="${dir}/../blastp_results"
      mkdir -p "$blast_output_dir"
      translated_file="${dir}/all_translated_sequences.faa.gz"
      output_blast="${blast_output_dir}/$(basename "${trimmed_file%.fastq.gz}_blastp.txt")"

      # Translate the sequences into proteins
      zcat "$trimmed_file" | seqtk seq -A - | seqkit translate -o - | awk 'NR%2==1 {header=$0} NR%2==0 {if (length($0) >= 30) print header; print $0}' | gzip > "$translated_file"

      # Run BLASTp
      zcat "$translated_file" | blastp -db "$PROTEIN_DB" -out "$output_blast" -outfmt 6 -evalue 1e-3 -num_threads 24
      echo "BLASTp output saved to $output_blast"

      # Extract matching sequence IDs and filter the FASTQ file
      awk '{print $1}' "$output_blast" | sort | uniq > "${blast_output_dir}/matched_ids.txt"
      seqtk subseq "$trimmed_file" "${blast_output_dir}/matched_ids.txt" > "${dir}/all_filtered_sequences.fastq.gz"
      echo "Filtered sequences saved to ${dir}/all_filtered_sequences.fastq.gz"
  done

  echo "BLASTp translation and check finished......."

  # Continue with further steps, such as alignment, only for filtered files
  echo "Alignment and BAM file generation in progress......."
  find "$FOLDER" -type f -path "*/all_filtered_sequences.fastq.gz" | while read -r filtered_file; do
      dir=$(dirname "$filtered_file")
      step1_dir="${dir}/../step1"
      mkdir -p "$step1_dir"
      basename=$(basename "$filtered_file" .fastq.gz)
      output_bam="$step1_dir/${basename}.bam"
      bwa mem -t 24 "$REFERENCE" "$filtered_file" | samtools view -Sb - > "$output_bam"
      echo "BAM file generated in step1 directory: $output_bam"
  done

  echo "Sorting and marking duplicates......."
  find "$FOLDER" -type f -path "*/step1/*.bam" | while read -r bam_file; do
      dir=$(dirname "$bam_file")
      basename=$(basename "$bam_file" .bam)
      step2_dir="${dir}/../step2"
      mkdir -p "$step2_dir"
      sorted_bam="${step2_dir}/${basename}_sorted.bam"
      marked_bam="${step2_dir}/${basename}_marked.bam"

      samtools sort -o "$sorted_bam" "$bam_file"
      echo "Sorted BAM file: $sorted_bam"

      samtools markdup "$sorted_bam" "$marked_bam"
      echo "Marked duplicates in BAM file: $marked_bam"
      # Index the BAM file to create BAI file
      samtools index "$sorted_bam"
      samtools index "$marked_bam"

      echo "BAI index file created for $output_bam"
  done

  echo "Sorting and marking duplicates finished......."

  echo "Feature counting with featureCounts......."
  find "$FOLDER" -mindepth 1 -maxdepth 1 -type d | while read -r r_folder; do
      r_group=$(basename "$r_folder")
      step2_dir="${r_folder}/step2"
      step3_dir="${r_folder}/step3"
      mkdir -p "$step3_dir"
      mapfile -t bam_files < <(find "$step2_dir" -type f -name "*.bam")

      if [ ${#bam_files[@]} -eq 0 ]; then
          echo "No BAM files found in $step2_dir, skipping."
          continue
      fi

      output_counts="${step3_dir}/${r_group}_combined_expression_counts.txt"
      featureCounts -a "$ANNOTATION" -o "$output_counts" -T 16 "${bam_files[@]}"
      echo "Combined feature counts for $r_group are in: $output_counts"
  done

  echo "All preprocessing done for folder $FOLDER. Please check the sequencing quality reports and existence of expression counts matrix file!"
done

