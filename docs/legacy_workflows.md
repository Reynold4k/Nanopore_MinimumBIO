# Legacy workflows and migration notes

The original pipelines are preserved verbatim under `legacy/` for provenance. This document explains how to run them as they were, maps each legacy file to its v2 replacement, and lists every intentional behaviour change between the legacy code and v2.

## Running the legacy pipelines as-is

The legacy scripts are **not configuration-driven**. To run them, edit the hard-coded paths inside each file to match your data and reference locations, then execute the relevant script directly:

```bash
cd legacy/Pipeline1
bash pipeline1.sh

# or for POD5 entry
bash pod5_preprocessing_with_genecounts.sh

# then run the R analysis after editing paths inside the R file
Rscript Pipeline1_limma.R
# or
Rscript Pipiline1_ttest.R
```

For the targeted workflow:

```bash
cd legacy/Pipeline1.5
bash pipeline1.5.sh
bash Compare.sh
Rscript Analysis.R
# or
Rscript Analysis_withoutidmapping.R
```

For candidate extraction and structural mapping:

```bash
cd legacy/Pipeline2
bash Pipeline2_1.sh      # genome mode
bash Pipeline2_1.5.sh    # targeted mode
bash Pipeline2_pymol.sh
```

## Legacy file → v2 script mapping

| Legacy file | v2 replacement | Notes |
|---|---|---|
| `Pipeline1/pipeline1.sh` | `workflow/preprocess_fastq.sh` | FASTQ genome-mode preprocessing |
| `Pipeline1/pod5_preprocessing_with_genecounts.sh` | `workflow/preprocess_pod5.sh` + `workflow/quantify_genome.sh` | POD5 basecalling + gene counts |
| `Pipeline1/Pipeline1_limma.R` | `R/differential_genome.R --method trend` | No limma modelling is actually performed (see below) |
| `Pipeline1/Pipiline1_ttest.R` | `R/differential_genome.R --method ttest` | Filename typo preserved in legacy |
| `Pipeline1.5/pipeline1.5.sh` + `Compare.sh` | `workflow/coverage_workflow.sh --stage preprocess/compare` | Targeted-reference coverage workflow |
| `Pipeline1.5/Analysis.R` | `R/differential_coverage.R --id-mapping <tsv>` | Adds UniProt id-mapping |
| `Pipeline1.5/Analysis_withoutidmapping.R` | `R/differential_coverage.R` | No id-mapping |
| `Pipeline2/Pipeline2_1.sh` | `workflow/candidate_extraction.sh` (genome mode) | Candidate sequence extraction |
| `Pipeline2/Pipeline2_1.5.sh` | `workflow/candidate_extraction.sh` (targeted mode) | Candidate extraction for targeted reference |
| `Pipeline2/Pipeline2_pymol.sh` | `workflow/structure_mapping.sh` | PDB alignment + PyMOL rendering |

## v2 behaviour changes vs legacy

The following changes are **documented, not silent**. Every item is traceable to `PROJECT_REVIEW.md` and `REFACTOR_PLAN.md` §4.

| Change | Legacy risk ID | New behaviour |
|---|---|---|
| C1 | R1 | `featureCounts` counts only the marked BAM; `--ignoreDup` added when markdup is enabled. |
| C2 | R2/R21 | Dorado is called as `dorado basecaller <model> <pod5_dir> > calls.bam`; `pod5 merge` is optional; `npmerge` is no longer used. |
| C3 | R3 | Startup check compares FASTA seqnames vs GTF seqnames and fails loudly on a chr-style mismatch. |
| C4 | R4 | Filtered FASTQ is genuinely gzipped. |
| C5 | R5 | The ≥30-aa length filter is actually applied (both `print` statements are inside the `awk` `if`). |
| C6 | R6 | Merge excludes `all_*` outputs; re-runs are safe with idempotent skip-if-exists logic. |
| C7 | R7 | Raw FASTQs are never deleted. |
| C8 | R8 | `differential_coverage.txt` is written with a header after sorting; R reads it with `header=TRUE`. |
| C9 | R9 | Experimental/control folders come only from config and are echoed at the start of every run. |
| C10 | R10 | Unified `pvalue_threshold` (default 0.05); a BH-adjusted p-value column is added. |
| C11 | R11/R12 | `bwa` uses `-x ont2d` (configurable); `featureCounts` uses `-L`. |
| C12 | R13 | Marking duplicates is optional via `alignment.mark_duplicates`. |
| C13 | R14 | `cutadapt` is used for motif filtering in targeted mode too, instead of the seqkit cross-line regex. |
| C14 | R15 | After BLASTp filtering, output is asserted to be non-empty; IDs are normalised with `awk '{print $1}'`. |
| C15 | R16 | GTF gene matching uses the `gene_name` attribute exactly, not whole-line substring. |
| C16 | R17 | Trend/t-test use all sample columns per round (mean over columns). |
| C17 | R18 | No fabricated rounds; trajectory plots require ≥2 rounds, otherwise they are skipped with a message. |
| C18 | R19/R20 | R scripts load only the packages they use: `edgeR`, `ggplot2`, `ggrepel`, `dplyr`, `tidyr`, `rtracklayer`, `readr`. |
| C19 | — | All scripts use `set -euo pipefail`, `require_tools`, and timestamped logging. |

### Preserved legacy behaviour

The following are intentionally kept as they were in the legacy scripts:

- Constant-region motif `GATCCGAATTCN`.
- BLASTp e-value default of `1e-3`.
- Top-1000 coverage loci in targeted mode.
- Growth-rate definition: last round minus first round, control-subtracted.
- PairwiseAligner scoring in structural alignment.
- PyMOL rendering style.
- Polynomial-trend and t-test methods themselves.

### Additional v2 deviations from legacy not in the table above

Two further v2 behaviours differ from the legacy scripts and are documented here:

1. **POD5 directory layout** (`workflow/preprocess_pod5.sh`):
   - Legacy `pod5_preprocessing_with_genecounts.sh` expected merged POD5 files and used `npmerge`.
   - v2 treats `input.pod5_dir` as a parent that contains round sub-directories. It accepts either `<pod5_dir>/<round>/` or `<pod5_dir>/pod5_pass/<round>/`, matching the standard Dorado pass/fail layout. Only pass reads are used.

2. **Candidate output layout** (`workflow/candidate_extraction.sh` and `workflow/structure_mapping.sh`):
   - Legacy scripts wrote candidate files directly under `<folder>/potential_hit/`.
   - v2 writes candidate and structure outputs under `<folder>/potential_hit/<round>/` so that results from multiple rounds do not collide.

## Important notes on legacy naming

- `Pipeline1_limma.R` is named after `limma`, but it does **not** perform limma modelling (`lmFit`/`eBayes`). It uses edgeR CPM, per-round differences, and polynomial trend regression. The v2 script is `R/differential_genome.R` with `--method trend`.
- `Pipiline1_ttest.R` contains a historical filename typo (`Pipiline` instead of `Pipeline`). It maps to `R/differential_genome.R --method ttest`.

## Recommendation

For new analyses, use the v2 workflow via `bin/minbio` with a YAML config file. Use the legacy scripts only when you need to reproduce an exact historical run or compare against the original implementation.
