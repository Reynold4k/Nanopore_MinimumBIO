#!/usr/bin/env bash
set -euo pipefail

# Genome-mode nanopore POD5 preprocessing.
# Modernised from legacy/Pipeline1/pod5_preprocessing_with_genecounts.sh.
#
# Dorado can basecall pass and fail reads separately. This workflow uses PASS
# reads only; point input.pod5_dir at the pass directory (or at a parent that
# contains round sub-directories of pass reads).
#
# Per experiment round: dorado basecaller -> BAM -> samtools fastq -> fastq.gz,
# then the work is delegated to workflow/preprocess_fastq.sh.

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
: "${MB_INPUT_POD5_DIR:?"missing input.pod5_dir"}"
: "${MB_INPUT_DORADO_MODEL:?"missing input.dorado_model"}"
: "${MB_INPUT_DORADO_DEVICE:?"missing input.dorado_device"}"

log "Experiment folder: $MB_EXPERIMENT_EXP_FOLDER"
log "Control folder:  $MB_EXPERIMENT_CONTROL_FOLDER"
log "Round pattern:     $MB_EXPERIMENT_ROUND_PATTERN"
log "POD5 dir:          $MB_INPUT_POD5_DIR"
log "Dorado model:      $MB_INPUT_DORADO_MODEL"
log "Dorado device:     $MB_INPUT_DORADO_DEVICE"

require_tools dorado samtools

# Resolve POD5 input for a round. Accepts either:
#   <pod5_dir>/<round>/
#   <pod5_dir>/pod5_pass/<round>/
find_pod5_dir() {
    local round_name="$1"
    local candidate
    for candidate in "$MB_INPUT_POD5_DIR/$round_name" "$MB_INPUT_POD5_DIR/pod5_pass/$round_name"; do
        if [[ -d "$candidate" ]]; then
            printf '%s' "$candidate"
            return 0
        fi
    done
    return 1
}

process_round() {
    local exp_round="$1"
    local round_name
    round_name=$(basename "$exp_round")
    local pod5_dir

    if ! pod5_dir=$(find_pod5_dir "$round_name"); then
        log "[$exp_round] No POD5 directory found for round $round_name; skipping"
        return 0
    fi

    local reports_dir="$exp_round/reports"
    mkdir -p "$reports_dir"

    local bam_file="$reports_dir/${round_name}_basecalled.bam"
    local fastq_file="$exp_round/all_sequences.fastq.gz"
    local flagstat_file="$reports_dir/${round_name}_flagstat.txt"

    if [[ -e "$bam_file" ]]; then
        log "[$exp_round] Basecalled BAM already exists: $bam_file (skipping)"
    else
        log "[$exp_round] Basecalling POD5: $pod5_dir"
        # C2: fixed dorado call emitting BAM to stdout.
        dorado basecaller \
            --device "$MB_INPUT_DORADO_DEVICE" \
            "$MB_INPUT_DORADO_MODEL" \
            "$pod5_dir" \
            > "$bam_file"
        log "[$exp_round] Basecalling complete: $bam_file"
    fi

    if [[ -e "$fastq_file" ]]; then
        log "[$exp_round] FASTQ already exists: $fastq_file (skipping)"
    else
        log "[$exp_round] Converting BAM to FASTQ..."
        samtools fastq "$bam_file" | gzip > "$fastq_file"
        log "[$exp_round] FASTQ complete: $fastq_file"
    fi

    if [[ -e "$flagstat_file" ]]; then
        log "[$exp_round] Flagstat report already exists: $flagstat_file (skipping)"
    else
        log "[$exp_round] Generating flagstat report..."
        samtools flagstat "$bam_file" > "$flagstat_file"
        log "[$exp_round] Flagstat complete: $flagstat_file"
    fi
}

if [[ ! -d "$MB_EXPERIMENT_EXP_FOLDER" ]]; then
    die "experiment folder does not exist: $MB_EXPERIMENT_EXP_FOLDER"
fi

while IFS= read -r -d '' exp_round; do
    process_round "$exp_round"
done < <(find "$MB_EXPERIMENT_EXP_FOLDER" -mindepth 1 -maxdepth 1 -type d -name "$MB_EXPERIMENT_ROUND_PATTERN" -print0)

log "Delegating to preprocess_fastq.sh..."
"$SCRIPT_DIR/preprocess_fastq.sh" --config "$config_file"

log "POD5 preprocessing complete."
