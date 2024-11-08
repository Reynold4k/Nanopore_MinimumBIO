#!/bin/bash


# This script automates the analysis and visualization of potential genomic hits from RNA sequencing data focused on specific genes.
# It is structured to process BAM files to extract regions of interest, convert these into FASTQ files, and perform a genome assembly
# using Canu for possible error correction. The process includes extracting hit regions from annotation files, converting data formats,
# and visualizing coverage. Overall, this pipeline is designed to facilitate gene-specific analyses within high-throughput datasets,
# allowing researchers to focus on particular genes of interest, manage raw data complexity, and enhance sequence accuracy with assembly strategies.



# What you need to modify:
# 1.Change the line 25 PARENT_DIR to your/path/to/exp/last_round/step2

# 2.Change the line 31 REFERENCE to your/reference/path

# 3.Change the line 32 ANNOTATION to your/ANNOTATION/path





# Paths to directories and files
# Base directory path; modify this line to change the directory structure
PARENT_DIR="/srv/scratch/z3546698/true/Small_Molecule/FK506/T7MB-2/231119/R3/step2"

# Derive the output and visualization directories from PARENT_DIR
OUTPUT_DIR="${PARENT_DIR%/*/*}/potential_hit"  # Remove the last two path segments (R3/step2)
VISUAL_DIR="${OUTPUT_DIR}/visualization"

ANN_FILE="/srv/scratch/z3546698/true/reference/Homo_sapiens.GRCh38.110.gtf"
GENOME_FASTA="/srv/scratch/z3546698/true/reference/hg38.fa"  # Reference FASTA file


FASTQ_DIR="${OUTPUT_DIR}/FASTQ"

# Ensure visualization directory exists
mkdir -p "$VISUAL_DIR"
mkdir -p "$FASTQ_DIR"
mkdir -p "$OUTPUT_DIR"

# Specify a list of target genes to analyze
GENES=("FKBP1A" "FKBP1C") # Gene names to search for and analyze within sequences

# Loop through each specified gene to extract data based on locations found in the GTF file
for gene_name in "${GENES[@]}"; do
  # Use awk to extract genomic locations for the gene's exons, prefixing with 'chr'
  Hit_LOCATIONS=$(awk -v gene="$gene_name" '$3 == "exon" && $0 ~ gene {print "chr"$1":"$4"-"$5}' "$ANN_FILE")

  # Skip if no location data is found for the specified gene in the annotation file
  if [ -z "$Hit_LOCATIONS" ]; then
    echo "No potential hit locations found for $gene_name in the annotation file."
    continue
  fi

  # Navigate to the directory containing trimmed and sorted BAM files
  cd "$PARENT_DIR"
  for bam_file in *_trimmed_sorted.bam; do
    # Ensure the BAM file has an index; create it if missing
    if [ ! -f "$bam_file.bai" ]; then
      echo "Indexing $bam_file..."
      samtools index "$bam_file"
    fi

    OUTPUT_BAM="${OUTPUT_DIR}/${gene_name}_Hit_${bam_file}"
    echo "Extracting potential hits for $gene_name from $bam_file to $OUTPUT_BAM..."

    # Extract reads for the specified locations using Samtools, saving to a new BAM file
    samtools view -b "$bam_file" $Hit_LOCATIONS > "$OUTPUT_BAM"
  done
done

echo "Extraction complete."

# Convert BAM files to FASTQ and run Canu individually per gene
cd "$OUTPUT_DIR"
for bam_file in *_Hit_*.bam; do
    # Extract the gene name from the BAM file name for further processing
    gene_name=$(basename "${bam_file}" | cut -d'_' -f1)
    
    # Convert the extracted BAM files containing potential hits to FASTQ format for Canu input
    fastq_file="${FASTQ_DIR}/${gene_name}.fastq"
    echo "Converting $bam_file to $fastq_file..."
    samtools fastq "$bam_file" > "$fastq_file"
    
    # Prepare an output directory for each gene's Canu assembly output
    canu_output_dir="${FASTQ_DIR}/canu_out_${gene_name}"
    mkdir -p "$canu_output_dir"
    
    # Run Canu for error correction or assembly of nanopore data
    echo "Running Canu assembly for $gene_name from $fastq_file..."
    canu \
        -p "${gene_name}_corrected" -d "$canu_output_dir" \
        genomeSize=0.038m \  # Estimated genome size for the assembly
        -nanopore-raw "$fastq_file" \  # Specify the FASTQ data format
        corOutCoverage=100           # Output coverage for correction
done

echo "Canu assembly for all genes complete."


# Use Bedtools to convert BAM to BED, sort, merge, and visualize coverage
cd "$OUTPUT_DIR"
for bam_file in *_Hit_*.bam; do

  # Convert BAM to BED format for downstream processing
  bed_file="${bam_file%.bam}.bed"
  bedtools bamtobed -i "$bam_file" > "$bed_file"
  
  # Sort the BED file by chromosome and start position
  sorted_bed_file="${bed_file%.bed}_sorted.bed"
  sort -k1,1 -k2,2n "$bed_file" > "$sorted_bed_file"

  # Merge overlapping regions and calculate coverage of merged intervals
  merged_bed_file="${bed_file%.bed}_merged.bed"
  bedtools merge -i "$sorted_bed_file" -c 4 -o count > "$merged_bed_file"

  # Convert BAM to BEDGraph for easier visualization of coverage metrics
  bedgraph_file="${bam_file%.bam}.bedgraph"
  bedtools genomecov -ibam "$bam_file" -bg > "$bedgraph_file"

  # If R is available, visualize coverage using ggplot2 and R
  if command -v Rscript &> /dev/null; then
    Rscript -e "
      library(ggplot2)
      library(dplyr)

      # Load BEDGraph data for coverage detail
      cov_data <- read.table('$bedgraph_file', header=FALSE)
      colnames(cov_data) <- c('chr', 'start', 'end', 'coverage')

      # Create a plot object with color gradient reflecting coverage
      p <- ggplot(cov_data, aes(x=start, xend=end, y=coverage, yend=coverage, fill=coverage)) +
        geom_segment(size=3) +  # Plot segments for coverage
        scale_fill_gradient(low='lightblue', high='darkred') +  # Define color gradient
        theme_minimal() +
        facet_wrap(~chr, scales='free_x', nrow=1) +  # Facet plot by chromosome
        labs(title='Coverage Visualization for ${bam_file}',
             x='Genomic Position',
             y='Coverage',
             fill='Coverage') +
        theme(strip.text.x = element_text(size = 8))

      # Save the plot to a file
      ggsave('${VISUAL_DIR}/${bam_file%.bam}_coverage_plot.png', plot=p, width=14, height=6)
    "
  else
    echo "Rscript is not available, skipping visualization for $bam_file"
  fi

done

echo "Analysis and visualization complete."

