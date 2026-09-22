#!/usr/bin/env bash
set -euo pipefail

# Gene-level quantification for genome mode.
# Runs featureCounts per round over the final sorted/marked BAM only (C1).

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

: "${MB_EXPERIMENT_EXP_FOLDER:?"missing experiment.exp_folder"}"
: "${MB_EXPERIMENT_CONTROL_FOLDER:?"missing experiment.control_folder"}"
: "${MB_EXPERIMENT_ROUND_PATTERN:?"missing experiment.round_pattern"}"
: "${MB_REFERENCE_ANNOTATION:?"missing reference.annotation"}"
: "${MB_ALIGNMENT_THREADS:?"missing alignment.threads"}"

mark_duplicates="${MB_ALIGNMENT_MARK_DUPLICATES:-true}"

log "Experiment folder: $MB_EXPERIMENT_EXP_FOLDER"
log "Control folder:  $MB_EXPERIMENT_CONTROL_FOLDER"
log "Annotation:        $MB_REFERENCE_ANNOTATION"
log "Threads:           $MB_ALIGNMENT_THREADS"
log "Mark duplicates:   $mark_duplicates"

require_tools featureCounts

quantify_round() {
    local round_dir="$1"
    local round_name
    round_name=$(basename "$round_dir")
    local step2_dir="$round_dir/step2"
    local step3_dir="$round_dir/step3"

    mkdir -p "$step3_dir"

    local output_counts="$step3_dir/${round_name}_combined_expression_counts.txt"
    if [[ -e "$output_counts" ]]; then
        log "[$round_dir] featureCounts output already exists: $output_counts (skipping)"
        return 0
    fi

    local bam_pattern
    if [[ "$mark_duplicates" == "true" ]]; then
        bam_pattern='*_marked.bam'
    else
        bam_pattern='*_sorted.bam'
    fi

    local bam_files=()
    while IFS= read -r -d '' bam; do
        bam_files+=("$bam")
    done < <(find "$step2_dir" -maxdepth 1 -type f -name "$bam_pattern" -print0)

    if [[ ${#bam_files[@]} -eq 0 ]]; then
        die "[$round_dir] no final BAM files matching $bam_pattern found in $step2_dir"
    fi

    local fc_opts=(-L -T "$MB_ALIGNMENT_THREADS" -a "$MB_REFERENCE_ANNOTATION" -o "$output_counts")
    if [[ "$mark_duplicates" == "true" ]]; then
        fc_opts+=(--ignoreDup)
    fi

    log "[$round_dir] Running featureCounts on ${#bam_files[@]} BAM(s)..."
    featureCounts "${fc_opts[@]}" "${bam_files[@]}"
    log "[$round_dir] featureCounts complete: $output_counts"
}

for base_folder in "$MB_EXPERIMENT_EXP_FOLDER" "$MB_EXPERIMENT_CONTROL_FOLDER"; do
    if [[ ! -d "$base_folder" ]]; then
        die "folder does not exist: $base_folder"
    fi

    log "Quantifying folder: $base_folder"
    while IFS= read -r -d '' round_dir; do
        quantify_round "$round_dir"
    done < <(find "$base_folder" -mindepth 1 -maxdepth 1 -type d -name "$MB_EXPERIMENT_ROUND_PATTERN" -print0)
done

log "Quantification complete."
