#!/usr/bin/env bash
set -euo pipefail

# Structural mapping: translated candidate proteins -> PDB alignment -> PyMOL.
# Modernised from legacy/Pipeline2/Pipeline2_pymol.sh.

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
: "${MB_STRUCTURE_GENES:?"missing structure.genes"}"
: "${MB_STRUCTURE_PDB_PATH:?"missing structure.pdb_path"}"

log "Experiment folder: $MB_EXPERIMENT_EXP_FOLDER"
log "Control folder:  $MB_EXPERIMENT_CONTROL_FOLDER"
log "PDB path:          $MB_STRUCTURE_PDB_PATH"
log "Genes:             $MB_STRUCTURE_GENES"

require_tools python3

pdb_align_py="$SCRIPT_DIR/../scripts/pdb_align.py"
if [[ ! -f "$pdb_align_py" ]]; then
    die "pdb_align.py not found: $pdb_align_py"
fi

pymol_available=false
if command -v pymol >/dev/null 2>&1; then
    pymol_available=true
fi

map_gene_to_protein_fa() {
    local folder="$1"
    local gene="$2"
    find "$folder/potential_hit" -type f -name "${gene}_translated_proteins.fa" -print0 2>/dev/null
}

process_protein() {
    local protein_fa="$1"
    local gene="$2"
    local out_dir="$(dirname "$protein_fa")/visualization"
    mkdir -p "$out_dir"

    local pml_file="$out_dir/${gene}_highlight.pml"
    local pse_file="$out_dir/${gene}_highlighted_structure.pse"
    local png_file="$out_dir/${gene}_highlighted_structure.png"

    if [[ -e "$png_file" ]]; then
        log "[$protein_fa] PNG already exists: $png_file (skipping)"
        return 0
    fi

    log "[$protein_fa] Aligning $gene against $MB_STRUCTURE_PDB_PATH..."
    local align_result
    if ! align_result=$(python3 "$pdb_align_py" --pdb "$MB_STRUCTURE_PDB_PATH" --query "$protein_fa" 2>/dev/null); then
        log "[$protein_fa] Alignment failed for $gene; skipping"
        return 0
    fi

    local chain_id start_index end_index
    read -r chain_id start_index end_index <<< "$align_result"

    if [[ -z "$chain_id" || -z "$start_index" || -z "$end_index" ]]; then
        log "[$protein_fa] Empty alignment result for $gene; skipping"
        return 0
    fi

    log "[$protein_fa] Best alignment: chain $chain_id, resi $start_index-$end_index"

    cat > "$pml_file" <<EOF
load $MB_STRUCTURE_PDB_PATH, structure
show cartoon, structure
select highlight_region, chain $chain_id and resi $start_index-$end_index
color gray, structure
color red, highlight_region
save $pse_file, format=pse
png $png_file, width=800, height=600, dpi=300, ray=1
quit
EOF
    log "[$protein_fa] PyMOL script: $pml_file"

    if [[ "$pymol_available" == "true" ]]; then
        log "[$protein_fa] Running PyMOL..."
        if pymol -c -d "$(cat "$pml_file")"; then
            log "[$protein_fa] PyMOL rendering complete: $png_file"
        else
            log "[$protein_fa] PyMOL failed; script saved at $pml_file"
            log "[$protein_fa] Run manually: pymol -c -q $pml_file"
        fi
    else
        log "[$protein_fa] PyMOL not available; script saved at $pml_file"
        log "[$protein_fa] Run manually: pymol -c -q $pml_file"
    fi
}

for folder in "$MB_EXPERIMENT_EXP_FOLDER" "$MB_EXPERIMENT_CONTROL_FOLDER"; do
    if [[ ! -d "$folder" ]]; then
        die "folder does not exist: $folder"
    fi

    log "Mapping structures for: $folder"
    for gene in $MB_STRUCTURE_GENES; do
        while IFS= read -r -d '' protein_fa; do
            process_protein "$protein_fa" "$gene"
        done < <(map_gene_to_protein_fa "$folder" "$gene")
    done
done

log "Structure mapping complete."
