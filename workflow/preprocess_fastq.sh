#!/usr/bin/env bash
set -euo pipefail

# Genome-mode nanopore FASTQ preprocessing.
# Modernised from legacy/Pipeline1/pipeline1.sh.
# Per round: merge -> cutadapt -> NanoPlot -> translate/BLASTp filter ->
# gzip -> bwa mem -> sort -> optional markdup -> index.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../scripts/lib/common.sh
source "$SCRIPT_DIR/../scripts/lib/common.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") --config <yaml>

Options:
  --config   Path to YAML configuration file (required).
EOF
}

config_file=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            if [[ $# -lt 2 ]]; then
                die "--config requires a value"
            fi
            config_file="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown argument: $1"
            ;;
    esac
done

if [[ -z "$config_file" ]]; then
    die "--config is required"
fi

load_config "$config_file"

# Validate required config values
: "${MB_EXPERIMENT_EXP_FOLDER:?"missing experiment.exp_folder"}"
: "${MB_EXPERIMENT_CONTROL_FOLDER:?"missing experiment.control_folder"}"
: "${MB_EXPERIMENT_ROUND_PATTERN:?"missing experiment.round_pattern"}"
: "${MB_REFERENCE_FASTA:?"missing reference.fasta"}"
: "${MB_REFERENCE_ANNOTATION:?"missing reference.annotation"}"
: "${MB_FILTERING_ADAPTER:?"missing filtering.adapter"}"
: "${MB_FILTERING_CUTADAPT_ERROR_RATE:?"missing filtering.cutadapt_error_rate"}"
: "${MB_FILTERING_MIN_PROTEIN_LENGTH:?"missing filtering.min_protein_length"}"
: "${MB_FILTERING_PROTEIN_DB:?"missing filtering.protein_db"}"
: "${MB_FILTERING_BLAST_EVALUE:?"missing filtering.blast_evalue"}"
: "${MB_ALIGNMENT_THREADS:?"missing alignment.threads"}"
: "${MB_ALIGNMENT_BWA_PRESET:?"missing alignment.bwa_preset"}"

# Boolean config values (default true for mark_duplicates)
mark_duplicates="${MB_ALIGNMENT_MARK_DUPLICATES:-true}"
remove_duplicates="${MB_ALIGNMENT_REMOVE_DUPLICATES:-false}"

log "Experiment folder: $MB_EXPERIMENT_EXP_FOLDER"
log "Control folder:  $MB_EXPERIMENT_CONTROL_FOLDER"
log "Round pattern:     $MB_EXPERIMENT_ROUND_PATTERN"
log "Reference:         $MB_REFERENCE_FASTA"
log "Annotation:        $MB_REFERENCE_ANNOTATION"
log "Protein DB:        $MB_FILTERING_PROTEIN_DB"
log "Threads:           $MB_ALIGNMENT_THREADS"
log "BWA preset:        $MB_ALIGNMENT_BWA_PRESET"
log "Mark duplicates:   $mark_duplicates"
log "Remove duplicates: $remove_duplicates"

require_tools cutadapt NanoPlot seqtk seqkit blastp gzip bwa samtools awk

process_round() {
    local round_dir="$1"
    local step1_dir="$round_dir/step1"
    local step2_dir="$round_dir/step2"
    local qc_dir="$round_dir/quality_control"
    local blast_dir="$round_dir/blastp_results"

    mkdir -p "$step1_dir" "$step2_dir" "$qc_dir" "$blast_dir"

    log "[$round_dir] Merging FASTQ files..."
    local merged_file="$round_dir/all_sequences.fastq.gz"
    if [[ -e "$merged_file" ]]; then
        log "[$round_dir] Merged file already exists: $merged_file (skipping)"
    else
        local temp_merged="${merged_file%.gz}"
        : > "$temp_merged"
        while IFS= read -r -d '' file; do
            log "[$round_dir] Merging file: $file"
            if [[ "$file" == *.gz ]]; then
                zcat "$file" >> "$temp_merged"
            else
                cat "$file" >> "$temp_merged"
            fi
        done < <(find "$round_dir" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) ! -name 'all_*' -print0)
        gzip "$temp_merged"
        log "[$round_dir] Merged and compressed: $merged_file"
    fi

    log "[$round_dir] Running cutadapt..."
    local trimmed_file="$step1_dir/all_trimmed.fastq.gz"
    if [[ -e "$trimmed_file" ]]; then
        log "[$round_dir] Trimmed file already exists: $trimmed_file (skipping)"
    else
        cutadapt -j "$MB_ALIGNMENT_THREADS" \
            -g "$MB_FILTERING_ADAPTER" \
            -e "$MB_FILTERING_CUTADAPT_ERROR_RATE" \
            --times=3 \
            --discard-untrimmed \
            -o "$trimmed_file" \
            "$merged_file"
        log "[$round_dir] Cutadapt complete: $trimmed_file"
    fi

    log "[$round_dir] Running NanoPlot..."
    local nanop_dir="$qc_dir/all_trimmed_nanop"
    if [[ -d "$nanop_dir" ]]; then
        log "[$round_dir] NanoPlot output already exists: $nanop_dir (skipping)"
    else
        NanoPlot --fastq "$trimmed_file" --outdir "$nanop_dir"
        log "[$round_dir] NanoPlot complete: $nanop_dir"
    fi

    log "[$round_dir] Translating and running BLASTp..."
    local translated_file="$step1_dir/all_translated_sequences.faa.gz"
    local output_blast="$blast_dir/all_trimmed_blastp.txt"
    local matched_ids="$blast_dir/matched_ids.txt"
    local filtered_file="$step1_dir/all_filtered_sequences.fastq.gz"

    if [[ -e "$filtered_file" ]]; then
        log "[$round_dir] Filtered file already exists: $filtered_file (skipping)"
    else
        if [[ -e "$translated_file" ]]; then
            log "[$round_dir] Translated file already exists: $translated_file (skipping)"
        else
            zcat "$trimmed_file" | seqtk seq -A - | seqkit translate -o - | \
                awk -v minlen="$MB_FILTERING_MIN_PROTEIN_LENGTH" \
                    'NR%2==1 {header=$0} NR%2==0 {if (length($0) >= minlen) {print header; print $0}}' | \
                gzip > "$translated_file"
            log "[$round_dir] Translation complete: $translated_file"
        fi

        if [[ -e "$output_blast" ]]; then
            log "[$round_dir] BLASTp output already exists: $output_blast (skipping)"
        else
            zcat "$translated_file" | blastp \
                -db "$MB_FILTERING_PROTEIN_DB" \
                -out "$output_blast" \
                -outfmt 6 \
                -evalue "$MB_FILTERING_BLAST_EVALUE" \
                -num_threads "$MB_ALIGNMENT_THREADS"
            log "[$round_dir] BLASTp complete: $output_blast"
        fi

        awk '{print $1}' "$output_blast" | sort -u > "$matched_ids"
        log "[$round_dir] Matched IDs: $matched_ids ($(wc -l < "$matched_ids") records)"

        seqtk subseq "$trimmed_file" "$matched_ids" | gzip > "$filtered_file"

        # C14: assert filtered FASTQ is non-empty
        if [[ ! -s "$filtered_file" ]]; then
            die "[$round_dir] filtered FASTQ is empty after BLASTp filtering: $filtered_file"
        fi
        log "[$round_dir] Filtered sequences complete: $filtered_file"
    fi

    log "[$round_dir] Aligning with bwa mem..."
    local basename="all_filtered_sequences"
    local raw_bam="$step1_dir/${basename}.bam"
    if [[ -e "$raw_bam" ]]; then
        log "[$round_dir] Raw BAM already exists: $raw_bam (skipping)"
    else
        bwa mem -x "$MB_ALIGNMENT_BWA_PRESET" -t "$MB_ALIGNMENT_THREADS" \
            "$MB_REFERENCE_FASTA" "$filtered_file" | \
            samtools view -Sb - > "$raw_bam"
        log "[$round_dir] Raw BAM complete: $raw_bam"
    fi

    log "[$round_dir] Sorting BAM..."
    local sorted_bam="$step2_dir/${basename}_sorted.bam"
    if [[ -e "$sorted_bam" ]]; then
        log "[$round_dir] Sorted BAM already exists: $sorted_bam (skipping)"
    else
        samtools sort -@ "$MB_ALIGNMENT_THREADS" -o "$sorted_bam" "$raw_bam"
        log "[$round_dir] Sorted BAM complete: $sorted_bam"
    fi

    local final_bam
    if [[ "$mark_duplicates" == "true" ]]; then
        log "[$round_dir] Marking duplicates..."
        local marked_bam="$step2_dir/${basename}_marked.bam"
        if [[ -e "$marked_bam" ]]; then
            log "[$round_dir] Marked BAM already exists: $marked_bam (skipping)"
            final_bam="$marked_bam"
        else
            local markdup_opts=()
            if [[ "$remove_duplicates" == "true" ]]; then
                markdup_opts+=(-r)
            fi
            samtools markdup "${markdup_opts[@]}" "$sorted_bam" "$marked_bam"
            log "[$round_dir] Marked BAM complete: $marked_bam"
            final_bam="$marked_bam"
        fi
    else
        final_bam="$sorted_bam"
    fi

    log "[$round_dir] Indexing final BAM: $final_bam"
    if [[ -e "${final_bam}.bai" ]]; then
        log "[$round_dir] BAI index already exists (skipping)"
    else
        samtools index "$final_bam"
        log "[$round_dir] Index complete: ${final_bam}.bai"
    fi
}

for base_folder in "$MB_EXPERIMENT_EXP_FOLDER" "$MB_EXPERIMENT_CONTROL_FOLDER"; do
    if [[ ! -d "$base_folder" ]]; then
        die "folder does not exist: $base_folder"
    fi

    log "Processing folder: $base_folder"
    while IFS= read -r -d '' round_dir; do
        process_round "$round_dir"
    done < <(find "$base_folder" -mindepth 1 -maxdepth 1 -type d -name "$MB_EXPERIMENT_ROUND_PATTERN" -print0)
done

log "Preprocessing complete."
