#!/usr/bin/env bash
set -euo pipefail

# Environment checker for Nanopore MinimumBIO v2.
# Usage: scripts/check_env.sh --config <config.yaml>

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/common.sh
source "$SCRIPT_DIR/lib/common.sh"

CONFIG_FILE=""

usage() {
    cat <<EOF
Usage: $(basename "$0") --config <config.yaml>

Checks that required tools, R packages and configured references are present.
Optional tools are reported but do not cause failure.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_FILE="${2:-}"
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

if [[ -z "$CONFIG_FILE" ]]; then
    usage >&2
    die "--config is required"
fi

log "Checking environment for MinimumBIO v2"

#------------------------------------------------------------------------------
# Required tools
#------------------------------------------------------------------------------
REQUIRED=(
    cutadapt
    seqkit
    seqtk
    blastp
    makeblastdb
    bwa
    minimap2
    samtools
    featureCounts
    NanoPlot
    bedtools
    python3
    Rscript
)

log "Checking required tools: ${REQUIRED[*]}"
MISSING=()
for tool in "${REQUIRED[@]}"; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        MISSING+=("$tool")
    fi
done
if ((${#MISSING[@]} > 0)); then
    log "MISSING required tools: ${MISSING[*]}"
fi

#------------------------------------------------------------------------------
# Optional tools (warn only)
#------------------------------------------------------------------------------
OPTIONAL=(dorado pymol pod5)
log "Checking optional tools: ${OPTIONAL[*]}"
for tool in "${OPTIONAL[@]}"; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        log "WARNING: optional tool not found: $tool"
    else
        log "OK: optional tool found: $tool"
    fi
done

#------------------------------------------------------------------------------
# R packages
#------------------------------------------------------------------------------
log "Checking R packages"
R_PACKAGES=(edgeR ggplot2 ggrepel dplyr tidyr rtracklayer readr)
if ! command -v Rscript >/dev/null 2>&1; then
    log "WARNING: Rscript not found; skipping R package check"
else
R_CHECK=$(cat <<EOF
missing <- c()
for (pkg in c("edgeR", "ggplot2", "ggrepel", "dplyr", "tidyr", "rtracklayer", "readr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) missing <- c(missing, pkg)
}
if (length(missing) > 0) {
  cat("MISSING_R:", paste(missing, collapse = " "), "\n", sep = "")
  quit(status = 1)
} else {
  cat("OK: all R packages present\n")
  quit(status = 0)
}
EOF
)

R_STATUS=0
R_OUTPUT=$(Rscript -e "$R_CHECK" 2>&1) || R_STATUS=$?
log "$R_OUTPUT"
if [[ "$R_STATUS" -ne 0 ]]; then
    MISSING+=("R packages")
fi
fi

#------------------------------------------------------------------------------
# Config-driven reference checks
#------------------------------------------------------------------------------
log "Loading config: $CONFIG_FILE"
load_config "$CONFIG_FILE"

FASTA="${MB_REFERENCE_FASTA:-}"
GTF="${MB_REFERENCE_ANNOTATION:-}"

if [[ -z "$FASTA" ]]; then
    die "config reference.fasta is not set"
fi

if [[ ! -f "$FASTA" ]]; then
    die "reference FASTA not found: $FASTA"
fi

log "Checking reference FASTA: $FASTA"

# Determine chr-prefix style of FASTA seqnames.
FASTA_PREFIX=none
if grep -qE '^>chr' "$FASTA" 2>/dev/null; then
    FASTA_PREFIX=chr
fi

if [[ -n "$GTF" && -f "$GTF" ]]; then
    log "Checking reference GTF: $GTF"

    # Extract first seqname from GTF (skip header lines).
    GTF_SEQNAME=$(awk '!/^#/ && NF>=1 {print $1; exit}' "$GTF")
    if [[ -z "$GTF_SEQNAME" ]]; then
        die "could not determine seqname from GTF: $GTF"
    fi

    GTF_PREFIX=none
    if [[ "$GTF_SEQNAME" == chr* ]]; then
        GTF_PREFIX=chr
    fi

    log "FASTA chr-prefix style: $FASTA_PREFIX"
    log "GTF chr-prefix style: $GTF_PREFIX (first seqname: $GTF_SEQNAME)"

    if [[ "$FASTA_PREFIX" != "$GTF_PREFIX" ]]; then
        die "chr-prefix mismatch between reference FASTA ($FASTA_PREFIX) and GTF ($GTF_PREFIX). Adjust reference files or use a consistent chr style."
    fi
else
    log "WARNING: no annotation GTF configured or file not found; skipping chr-prefix check"
fi

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------
if ((${#MISSING[@]} > 0)); then
    die "environment check failed. Missing: ${MISSING[*]}"
fi

log "Environment check passed"
