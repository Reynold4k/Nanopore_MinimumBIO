#!/usr/bin/env bash
set -euo pipefail

# Candidate sequence extraction: BAM/GTF -> BED -> DNA -> protein.
# Merges legacy/Pipeline2/Pipeline2_1.sh and Pipeline2_1.5.sh.

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
: "${MB_MODE:?"missing mode"}"
: "${MB_REFERENCE_FASTA:?"missing reference.fasta"}"
: "${MB_STRUCTURE_GENES:?"missing structure.genes"}"
: "${MB_STRUCTURE_CHR_PREFIX:-auto}"

log "Experiment folder: $MB_EXPERIMENT_EXP_FOLDER"
log "Control folder:  $MB_EXPERIMENT_CONTROL_FOLDER"
log "Mode:              $MB_MODE"
log "Reference:         $MB_REFERENCE_FASTA"
log "Genes:             $MB_STRUCTURE_GENES"
log "Chr prefix:        ${MB_STRUCTURE_CHR_PREFIX:-auto}"

require_tools samtools bedtools python3 awk

translate_dna_py="$SCRIPT_DIR/../scripts/translate_dna.py"
if [[ ! -f "$translate_dna_py" ]]; then
    die "translate_dna.py not found: $translate_dna_py"
fi

# Detect chromosome naming style from the first @SQ line of a BAM.
detect_chr_prefix() {
    local bam="$1"
    local sq
    sq=$(samtools view -H "$bam" | grep '^@SQ' | head -n 1)
    if [[ "$sq" =~ SN:chr ]]; then
        printf 'chr'
    else
        printf 'none'
    fi
}

# Extract gene_name value from a GTF attributes field and compare to target.
extract_exons_genome() {
    local gtf="$1"
    local gene="$2"
    awk -v gene="$gene" -F'\t' '
        $3 == "exon" {
            attrs = $9
            # Strip leading/trailing spaces and quotes around values.
            gsub(/^\s+|\s+$/, "", attrs)
            split(attrs, parts, ";")
            for (i in parts) {
                gsub(/^\s+|\s+$/, "", parts[i])
                if (parts[i] ~ /^gene_name\s+/) {
                    val = parts[i]
                    sub(/^gene_name\s+/, "", val)
                    gsub(/^"|"$/, "", val)
                    if (val == gene) {
                        print $1 "\t" $4 "\t" $5
                        break
                    }
                }
            }
        }
    ' "$gtf"
}

extract_genome_round() {
    local folder="$1"
    local round_dir="$2"
    local round_name
    round_name=$(basename "$round_dir")
    local step2_dir="$round_dir/step2"
    local out_dir="$folder/potential_hit/$round_name"
    mkdir -p "$out_dir"

    local bam_files=()
    while IFS= read -r -d '' bam; do
        bam_files+=("$bam")
    done < <(find "$step2_dir" -maxdepth 1 -type f -name '*_sorted.bam' -print0)

    if [[ ${#bam_files[@]} -eq 0 ]]; then
        log "[$round_dir] No sorted BAM found; skipping"
        return 0
    fi

    local bam="${bam_files[0]}"
    log "[$round_dir] Using BAM: $bam"

    local chr_prefix="${MB_STRUCTURE_CHR_PREFIX:-auto}"
    if [[ "$chr_prefix" == "auto" ]]; then
        chr_prefix=$(detect_chr_prefix "$bam")
        log "[$round_dir] Auto-detected chr prefix: $chr_prefix"
    fi

    for gene in $MB_STRUCTURE_GENES; do
        log "[$round_dir] Extracting candidate regions for gene: $gene"

        local exons_bed="$out_dir/${gene}_exons.bed"
        if [[ -e "$exons_bed" ]]; then
            log "[$round_dir] Exons BED already exists: $exons_bed (skipping)"
        else
            extract_exons_genome "$MB_REFERENCE_ANNOTATION" "$gene" > "$exons_bed"
            log "[$round_dir] Exons BED: $exons_bed"
        fi

        if [[ ! -s "$exons_bed" ]]; then
            log "[$round_dir] No exons found for $gene; skipping"
            continue
        fi

        local regions=()
        while read -r seqname start end; do
            local region_seqname="$seqname"
            if [[ "$chr_prefix" == "chr" && "$seqname" != chr* ]]; then
                region_seqname="chr${seqname}"
            elif [[ "$chr_prefix" == "none" ]]; then
                region_seqname="${seqname#chr}"
            fi
            regions+=("${region_seqname}:${start}-${end}")
        done < "$exons_bed"

        local hit_bam="$out_dir/${gene}_Hit_$(basename "$bam")"
        if [[ -e "$hit_bam" ]]; then
            log "[$round_dir] Hit BAM already exists: $hit_bam (skipping)"
        else
            samtools view -b "$bam" "${regions[@]}" > "$hit_bam"
            log "[$round_dir] Hit BAM: $hit_bam"
        fi

        local bed_file="${hit_bam%.bam}.bed"
        local sorted_bed="${bed_file%.bed}_sorted.bed"
        if [[ -e "$sorted_bed" ]]; then
            log "[$round_dir] Sorted BED already exists: $sorted_bed (skipping)"
        else
            bedtools bamtobed -i "$hit_bam" > "$bed_file"
            sort -k1,1 -k2,2n "$bed_file" > "$sorted_bed"
            log "[$round_dir] Sorted BED: $sorted_bed"
        fi

        local dna_out="$out_dir/${gene}_high_coverage_sequences.fa"
        local protein_out="$out_dir/${gene}_translated_proteins.fa"
        if [[ -e "$protein_out" ]]; then
            log "[$round_dir] Translated proteins already exist: $protein_out (skipping)"
            continue
        fi

        bedtools getfasta -fi "$MB_REFERENCE_FASTA" -bed "$sorted_bed" -fo "$dna_out"
        log "[$round_dir] DNA sequences: $dna_out"

        python3 "$translate_dna_py" -i "$dna_out" -o "$protein_out"
        log "[$round_dir] Translated proteins: $protein_out"
    done
}

extract_targeted_round() {
    local folder="$1"
    local round_dir="$2"
    local vis_dir="$round_dir/step2/visualization"
    local out_dir="$folder/potential_hit/$(basename "$round_dir")"
    mkdir -p "$out_dir"

    if [[ ! -d "$vis_dir" ]]; then
        log "[$round_dir] No visualization BED directory: $vis_dir; skipping"
        return 0
    fi

    local bed_count=0
    while IFS= read -r -d '' bed_file; do
        bed_count=$((bed_count + 1))
        local base
        base=$(basename "$bed_file" .bed)
        local dna_out="$out_dir/${base}_high_coverage_sequences.fa"
        local protein_out="$out_dir/${base}_translated_proteins.fa"

        if [[ -e "$protein_out" ]]; then
            log "[$round_dir] Translated proteins already exist: $protein_out (skipping)"
            continue
        fi

        bedtools getfasta -fi "$MB_REFERENCE_FASTA" -bed "$bed_file" -fo "$dna_out"
        log "[$round_dir] DNA sequences for $base: $dna_out"

        python3 "$translate_dna_py" -i "$dna_out" -o "$protein_out"
        log "[$round_dir] Translated proteins for $base: $protein_out"
    done < <(find "$vis_dir" -maxdepth 1 -type f -name '*.bed' -print0)

    if [[ "$bed_count" -eq 0 ]]; then
        log "[$round_dir] No BED files found in $vis_dir"
    fi
}

for folder in "$MB_EXPERIMENT_EXP_FOLDER" "$MB_EXPERIMENT_CONTROL_FOLDER"; do
    if [[ ! -d "$folder" ]]; then
        die "folder does not exist: $folder"
    fi

    log "Extracting candidates from: $folder"
    while IFS= read -r -d '' round_dir; do
        case "$MB_MODE" in
            genome)
                : "${MB_REFERENCE_ANNOTATION:?"missing reference.annotation for genome mode"}"
                extract_genome_round "$folder" "$round_dir"
                ;;
            targeted)
                extract_targeted_round "$folder" "$round_dir"
                ;;
            *)
                die "unknown mode: $MB_MODE"
                ;;
        esac
    done < <(find "$folder" -mindepth 1 -maxdepth 1 -type d -name "$MB_EXPERIMENT_ROUND_PATTERN" -print0)
done

log "Candidate extraction complete."
