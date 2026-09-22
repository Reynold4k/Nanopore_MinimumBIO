#!/usr/bin/env bash
set -euo pipefail

# Shared shell library for Nanopore MinimumBIO v2.
# Source this file in other scripts: source "$(dirname "$0")/lib/common.sh"

#------------------------------------------------------------------------------
# Logging helpers
#------------------------------------------------------------------------------

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

die() {
    log "ERROR: $*"
    exit 1
}

#------------------------------------------------------------------------------
# Tool dependency check
#------------------------------------------------------------------------------

require_tools() {
    local missing=()
    local tool
    for tool in "$@"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            missing+=("$tool")
        fi
    done
    if ((${#missing[@]} > 0)); then
        die "missing required tools: ${missing[*]}"
    fi
}

#------------------------------------------------------------------------------
# YAML config loader
# Exports flat shell variables named MB_<SECTION>_<KEY> (uppercased).
# Lists are exported as space-separated strings.
#------------------------------------------------------------------------------

load_config() {
    local config_file="${1:-}"
    if [[ -z "$config_file" ]]; then
        die "load_config: no config file provided"
    fi
    if [[ ! -f "$config_file" ]]; then
        die "load_config: config file not found: $config_file"
    fi
    if ! command -v python3 >/dev/null 2>&1; then
        die "load_config: python3 is required to parse YAML config"
    fi
    if ! python3 -c 'import yaml' >/dev/null 2>&1; then
        die "load_config: PyYAML is required to parse YAML config (install with: conda install pyyaml)"
    fi

    local python_script
    python_script=$(cat <<'PY'
import sys, yaml

def export_flat(prefix, obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            export_flat(f"{prefix}_{k.upper()}", v)
    elif isinstance(obj, list):
        # Space-separated string; quote each element to protect whitespace.
        vals = ' '.join(str(x) for x in obj)
        print(f"export {prefix}={vals!r}")
    elif isinstance(obj, bool):
        print(f"export {prefix}={'true' if obj else 'false'}")
    elif obj is None:
        print(f'export {prefix}=""')
    else:
        print(f"export {prefix}={str(obj)!r}")

with open(sys.argv[1], 'r') as fh:
    cfg = yaml.safe_load(fh)
if not isinstance(cfg, dict):
    print(f"Top-level config is not a mapping", file=sys.stderr)
    sys.exit(1)
for section, value in cfg.items():
    export_flat(f"MB_{section.upper()}", value)
PY
)

    local exports
    if ! exports=$(python3 -c "$python_script" "$config_file" 2>&1); then
        die "load_config: failed to parse YAML config: $exports"
    fi
    # shellcheck disable=SC1090
    eval "$exports"
}
