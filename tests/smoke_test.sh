#!/usr/bin/env bash
set -euo pipefail

# Smoke test for Nanopore MinimumBIO v2.
# Validates shell/R syntax and the example config without requiring the heavy
# bioinformatics tools to be installed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

FAILED=0
PASS=0

pass() {
    echo "  PASS: $*" >&2
    ((PASS++)) || true
}

fail() {
    echo "  FAIL: $*" >&2
    ((FAILED++)) || true
}

#------------------------------------------------------------------------------
# 1. Shell script syntax
#------------------------------------------------------------------------------
echo "Checking shell script syntax ..." >&2

for dir in workflow scripts bin; do
    dir_path="$REPO_ROOT/$dir"
    if [[ ! -d "$dir_path" ]]; then
        echo "  SKIP: directory does not exist: $dir" >&2
        continue
    fi
    while IFS= read -r -d '' script; do
        if bash -n "$script"; then
            pass "bash -n $script"
        else
            fail "bash -n $script"
        fi
    done < <(find "$dir_path" -type f -name '*.sh' -print0)
done

#------------------------------------------------------------------------------
# 2. R script syntax
#------------------------------------------------------------------------------
echo "Checking R script syntax ..." >&2

if command -v Rscript >/dev/null 2>&1; then
    for dir in R legacy; do
        dir_path="$REPO_ROOT/$dir"
        if [[ ! -d "$dir_path" ]]; then
            echo "  SKIP: directory does not exist: $dir" >&2
            continue
        fi
        while IFS= read -r -d '' rscript; do
            if Rscript -e "invisible(parse(\"$rscript\"))" >/dev/null 2>&1; then
                pass "R parse $rscript"
            else
                fail "R parse $rscript"
            fi
        done < <(find "$dir_path" -type f \( -name '*.R' -o -name '*.r' \) -print0)
    done
else
    echo "  SKIP: Rscript not available" >&2
fi

#------------------------------------------------------------------------------
# 3. YAML config validation (same logic as scripts/lib/common.sh)
#------------------------------------------------------------------------------
echo "Checking YAML config parsing ..." >&2

if command -v python3 >/dev/null 2>&1 && python3 -c 'import yaml' >/dev/null 2>&1; then
    python_script=$(cat <<'PY'
import sys, yaml

def export_flat(prefix, obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            export_flat(f"{prefix}_{k.upper()}", v)
    elif isinstance(obj, list):
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
    print("Top-level config is not a mapping", file=sys.stderr)
    sys.exit(1)
for section, value in cfg.items():
    export_flat(f"MB_{section.upper()}", value)
PY
)

    for cfg in "$REPO_ROOT/config/config.example.yaml" "$REPO_ROOT/examples/config.fastq.yaml" "$REPO_ROOT/examples/config.pod5.yaml" "$REPO_ROOT/examples/config.targeted.yaml"; do
        if [[ -f "$cfg" ]]; then
            if python3 -c "$python_script" "$cfg" >/dev/null 2>&1; then
                pass "YAML parse $cfg"
            else
                fail "YAML parse $cfg"
            fi
        else
            fail "missing config: $cfg"
        fi
    done
else
    echo "  SKIP: python3 or PyYAML not available" >&2
fi

#------------------------------------------------------------------------------
# 4. Absolute path check
#------------------------------------------------------------------------------
echo "Checking for hard-coded absolute paths ..." >&2

# Find files under repo root that are not in legacy/ and not the planning docs.
while IFS= read -r -d '' file; do
    case "$file" in
        */legacy/*|*/PROJECT_REVIEW.md|*/REFACTOR_PLAN.md|*/docs/*)
            continue
            ;;
    esac
    if grep -nE '(^|[[:space:]])(/mnt/|/srv/scratch|D:/)' "$file" >/dev/null 2>&1; then
        fail "absolute path in $file"
    else
        pass "no absolute paths in $file"
    fi
done < <(find "$REPO_ROOT" -type f \( -name '*.sh' -o -name '*.R' -o -name '*.r' -o -name '*.yaml' -o -name '*.yml' -o -name '*.md' -o -name '*.py' \) -print0)

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------
echo "" >&2
echo "Smoke test summary: $PASS passed, $FAILED failed" >&2
if [[ "$FAILED" -gt 0 ]]; then
    exit 1
fi
exit 0
