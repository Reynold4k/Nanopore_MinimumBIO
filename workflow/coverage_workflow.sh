#!/usr/bin/env bash
set -euo pipefail

# Targeted-reference coverage workflow.
# Modernised from legacy/Pipeline1.5/pipeline1.5.sh + Compare.sh.
#
# Usage:
#   coverage_workflow.sh --config <yaml> --stage preprocess
#   coverage_workflow.sh --config <yaml> --stage compare

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../scripts/lib/common.sh
source "$SCRIPT_DIR/../scripts/lib/common.sh"

usage() {
    cat <<EOF
Usage: $(basename "$0") --config <yaml> --stage <preprocess|compare>

Options:
  --config   Path to YAML configuration file (required).
  --stage    Stage to run: preprocess or compare (required).
EOF
}

config_file=""
stage=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            if [[ $# -lt 2 ]]; then
                die "--config requires a value"
            fi
            config_file="$2"
            shift 2
            ;;
        --stage)
            if [[ $# -lt 2 ]]; then
                die "--stage requires a value"
            fi
            stage="$2"
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
if [[ -z "$stage" ]]; then
    die "--stage is required"
fi
if [[ "$stage" != "preprocess" && "$stage" != "compare" ]]; then
    die "--stage must be 'preprocess' or 'compare'"
fi

load_config "$config_file"

: "${MB_EXPERIMENT_EXP_FOLDER:?"missing experiment.exp_folder"}"
: "${MB_EXPERIMENT_CONTROL_FOLDER:?"missing experiment.control_folder"}"
: "${MB_EXPERIMENT_ROUND_PATTERN:?"missing experiment.round_pattern"}"
: "${MB_REFERENCE_FASTA:?"missing reference.fasta"}"
: "${MB_FILTERING_ADAPTER:?"missing filtering.adapter"}"
: "${MB_FILTERING_CUTADAPT_ERROR_RATE:?"missing filtering.cutadapt_error_rate"}"
: "${MB_FILTERING_MIN_PROTEIN_LENGTH:?"missing filtering.min_protein_length"}"
: "${MB_FILTERING_PROTEIN_DB:?"missing filtering.protein_db"}"
: "${MB_FILTERING_BLAST_EVALUE:?"missing filtering.blast_evalue"}"
: "${MB_ALIGNMENT_THREADS:?"missing alignment.threads"}"

log "Experiment folder: $MB_EXPERIMENT_EXP_FOLDER"
log "Control folder:  $MB_EXPERIMENT_CONTROL_FOLDER"
log "Round pattern:     $MB_EXPERIMENT_ROUND_PATTERN"
log "Reference:         $MB_REFERENCE_FASTA"
log "Stage:             $stage"

require_tools cutadapt NanoPlot seqtk seqkit blastp gzip minimap2 samtools awk

merge_round() {
    local round_dir="$1"
    local merged_file="$round_dir/all_sequences.fastq.gz"

    if [[ -e "$merged_file" ]]; then
        log "[$round_dir] Merged file already exists: $merged_file (skipping)"
        return 0
    fi

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
}

filter_and_align_round() {
    local round_dir="$1"
    local step1_dir="$round_dir/step1"
    local step2_dir="$round_dir/step2"
    local qc_dir="$round_dir/quality_control"
    local blast_dir="$round_dir/blastp_results"

    mkdir -p "$step1_dir" "$step2_dir" "$qc_dir" "$blast_dir"

    merge_round "$round_dir"

    local merged_file="$round_dir/all_sequences.fastq.gz"
    local trimmed_file="$step1_dir/all_trimmed.fastq.gz"

    if [[ -e "$trimmed_file" ]]; then
        log "[$round_dir] Trimmed file already exists: $trimmed_file (skipping)"
    else
        log "[$round_dir] Running cutadapt motif filter..."
        cutadapt -j "$MB_ALIGNMENT_THREADS" \
            -g "$MB_FILTERING_ADAPTER" \
            -e "$MB_FILTERING_CUTADAPT_ERROR_RATE" \
            --times=3 \
            --discard-untrimmed \
            -o "$trimmed_file" \
            "$merged_file"
        log "[$round_dir] Cutadapt complete: $trimmed_file"
    fi

    local nanop_dir="$qc_dir/all_trimmed_nanop"
    if [[ -d "$nanop_dir" ]]; then
        log "[$round_dir] NanoPlot output already exists: $nanop_dir (skipping)"
    else
        log "[$round_dir] Running NanoPlot..."
        NanoPlot --fastq "$trimmed_file" --outdir "$nanop_dir"
        log "[$round_dir] NanoPlot complete: $nanop_dir"
    fi

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
        seqtk subseq "$trimmed_file" "$matched_ids" | gzip > "$filtered_file"

        if [[ ! -s "$filtered_file" ]]; then
            die "[$round_dir] filtered FASTQ is empty after BLASTp filtering: $filtered_file"
        fi
        log "[$round_dir] Filtered sequences complete: $filtered_file"
    fi

    local bam_file="$step2_dir/aligned.bam"
    if [[ -e "$bam_file" ]]; then
        log "[$round_dir] Aligned BAM already exists: $bam_file (skipping)"
    else
        log "[$round_dir] Aligning with minimap2..."
        minimap2 -ax map-ont -t "$MB_ALIGNMENT_THREADS" "$MB_REFERENCE_FASTA" "$filtered_file" | \
            samtools view -Sb - > "$bam_file"
        log "[$round_dir] Alignment complete: $bam_file"
    fi

    local sorted_bam="$step2_dir/aligned_sorted.bam"
    if [[ -e "$sorted_bam" ]]; then
        log "[$round_dir] Sorted BAM already exists: $sorted_bam (skipping)"
    else
        log "[$round_dir] Sorting BAM..."
        samtools sort -@ "$MB_ALIGNMENT_THREADS" -o "$sorted_bam" "$bam_file"
        log "[$round_dir] Sort complete: $sorted_bam"
    fi

    if [[ -e "${sorted_bam}.bai" ]]; then
        log "[$round_dir] BAI index already exists (skipping)"
    else
        samtools index "$sorted_bam"
        log "[$round_dir] Index complete: ${sorted_bam}.bai"
    fi

    local coverage_file="$step2_dir/coverage.txt"
    if [[ -e "$coverage_file" ]]; then
        log "[$round_dir] Coverage file already exists: $coverage_file (skipping)"
    else
        log "[$round_dir] Computing coverage..."
        samtools depth "$sorted_bam" > "$coverage_file"
        log "[$round_dir] Coverage complete: $coverage_file"
    fi

    local positions_file="$step2_dir/positions.txt"
    local top_file="$step2_dir/top_1000_positions.txt"
    if [[ -e "$top_file" ]]; then
        log "[$round_dir] Top-1000 positions already exist: $top_file (skipping)"
        return 0
    fi

    log "[$round_dir] Aggregating coverage by reference..."
    awk '{
        key = $1
        position = $2
        coverage = $3
        if (key in positions) {
            if (position < positions[key][1]) {
                positions[key][1] = position
            }
            positions[key][2] = position
            positions[key][3] += coverage
        } else {
            positions[key][1] = position
            positions[key][2] = position
            positions[key][3] = coverage
        }
    } END {
        for (gene in positions) {
            print gene "\t" positions[gene][1] "\t" positions[gene][2] "\t" positions[gene][3]
        }
    }' "$coverage_file" > "$positions_file"

    sort -k4,4nr "$positions_file" | head -n 1000 > "$top_file"
    log "[$round_dir] Top-1000 positions: $top_file"

    local vis_dir="$step2_dir/visualization"
    mkdir -p "$vis_dir"

    local counter=0
    local max_bed=1000
    while IFS=$'\t' read -r gene_id start end total_coverage; do
        if [[ "$counter" -ge "$max_bed" ]]; then
            break
        fi
        local bed_file="$vis_dir/${gene_id}.bed"
        if [[ -e "$bed_file" ]]; then
            log "[$round_dir] BED already exists: $bed_file (skipping)"
        else
            printf '%s\t%s\t%s\n' "$gene_id" "$start" "$end" > "$bed_file"
            log "[$round_dir] Created BED: $bed_file"
        fi
        counter=$((counter + 1))
    done < "$top_file"
}

compare_rounds() {
    if [[ ! -d "$MB_EXPERIMENT_EXP_FOLDER" ]]; then
        die "experiment folder does not exist: $MB_EXPERIMENT_EXP_FOLDER"
    fi
    if [[ ! -d "$MB_EXPERIMENT_CONTROL_FOLDER" ]]; then
        die "control folder does not exist: $MB_EXPERIMENT_CONTROL_FOLDER"
    fi

    while IFS= read -r -d '' exp_round; do
        local round_name
        round_name=$(basename "$exp_round")
        local ctrl_round="$MB_EXPERIMENT_CONTROL_FOLDER/$round_name"

        if [[ ! -d "$ctrl_round" ]]; then
            log "[$exp_round] No matching control round $ctrl_round; skipping"
            continue
        fi

        local top_exp="$exp_round/step2/top_1000_positions.txt"
        local top_ctrl="$ctrl_round/step2/top_1000_positions.txt"

        if [[ ! -f "$top_exp" || ! -f "$top_ctrl" ]]; then
            log "[$exp_round] Missing top-1000 file(s); skipping"
            continue
        fi

        local output_file="$exp_round/differential_coverage.txt"
        if [[ -e "$output_file" ]]; then
            log "[$exp_round] Differential coverage already exists: $output_file (skipping)"
            continue
        fi

        log "[$exp_round] Comparing coverage with $ctrl_round..."

        declare -A control_coverage
        declare -A experimental_coverage

        while read -r gene_id _start _end coverage; do
            control_coverage["$gene_id"]="$coverage"
        done < <(cut -f 1,2,3,4 "$top_ctrl")

        while read -r gene_id _start _end coverage; do
            experimental_coverage["$gene_id"]=$(awk "BEGIN {print (${experimental_coverage[$gene_id]:-0} + $coverage)}")
        done < <(cut -f 1,2,3,4 "$top_exp")

        local tmp_body
        tmp_body=$(mktemp)
        # shellcheck disable=SC2064
        trap "rm -f '$tmp_body'" RETURN

        for gene_id in "${!experimental_coverage[@]}"; do
            local ctrl_cov="${control_coverage[$gene_id]:-0}"
            local exp_cov="${experimental_coverage[$gene_id]:-0}"
            local difference
            difference=$(awk "BEGIN {print $exp_cov - $ctrl_cov}")
            if awk "BEGIN {exit !($difference != 0)}"; then
                printf '%s\t%s\t%s\t%s\n' "$gene_id" "$ctrl_cov" "$exp_cov" "$difference" >> "$tmp_body"
            fi
        done

        # C8: write header after sorting so it stays on top; use temp file, not sort -o.
        sort -k4,4gr "$tmp_body" > "$output_file"
        {
            printf 'Gene\tControl_Coverage\tExperimental_Coverage\tDifference\n'
            cat "$output_file"
        } > "${output_file}.tmp"
        mv "${output_file}.tmp" "$output_file"

        log "[$exp_round] Differential coverage complete: $output_file"
    done < <(find "$MB_EXPERIMENT_EXP_FOLDER" -mindepth 1 -maxdepth 1 -type d -name "$MB_EXPERIMENT_ROUND_PATTERN" -print0)
}

case "$stage" in
    preprocess)
        if [[ ! -d "$MB_EXPERIMENT_EXP_FOLDER" ]]; then
            die "experiment folder does not exist: $MB_EXPERIMENT_EXP_FOLDER"
        fi
        if [[ ! -d "$MB_EXPERIMENT_CONTROL_FOLDER" ]]; then
            die "control folder does not exist: $MB_EXPERIMENT_CONTROL_FOLDER"
        fi
        for base_folder in "$MB_EXPERIMENT_EXP_FOLDER" "$MB_EXPERIMENT_CONTROL_FOLDER"; do
            log "Preprocessing folder: $base_folder"
            while IFS= read -r -d '' round_dir; do
                filter_and_align_round "$round_dir"
            done < <(find "$base_folder" -mindepth 1 -maxdepth 1 -type d -name "$MB_EXPERIMENT_ROUND_PATTERN" -print0)
        done
        log "Preprocessing stage complete."
        ;;
    compare)
        compare_rounds
        log "Compare stage complete."
        ;;
esac
