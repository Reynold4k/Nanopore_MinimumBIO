# Project history

This narrative is reconstructed from the repository's git history and from `PROJECT_REVIEW.md` §5. It frames the codebase as experiment-driven iterative research software.

## Timeline

### 2024-09 — Initial upload and HPC adaptation

- **17 Sep 2024**: initial upload of `porechop_preprocessing.sh`, Dorado usage, `Pipeline1.R`, and the first version of Pipeline2 (`Pipeline2_a.sh` / `Pipeline2_b.sh`).
- **18 Sep 2024**: HPC adaptation for the Katana cluster, including `module load` versions.
- **19 Sep 2024**: Pipeline2 gained BLAST and AlphaFold/PyMOL ideas.
- **22 Sep 2024**: added `pod5_preprocessing_with_genecounts.sh` for POD5 entry.
- **27 Sep 2024**: **Pipeline1.5 born** — the targeted-reference coverage workflow, plus `Compare.sh` for experimental/control comparison.
- **30 Sep 2024**: `Analysis.R` and `Pipeline2_1.5.sh` added, giving the personalised/targeted branch of Pipeline2.

September was the most active month, with **232 commits** as the core workflow took shape.

### 2024-10 — Preprocessing method iterations

- **9 Oct 2024**: `seqkit` motif regex introduced (`GATCCGAATTC[ACGTN]`) for adapter/constant-region filtering.
- **30 Oct 2024**: Porechop dropped from the main FASTQ path.

October had **196 commits**. The preprocessing approach was still being refined.

### 2024-11 — Protein-aware filtering and analysis split

- **4 Nov 2024**: BLASTx in-frame check added.
- **7–8 Nov 2024**: `.pbs` scripts migrated to `.sh` for portability off the HPC.
- **15 Nov 2024**: Pipeline2 consolidated into `Pipeline2_1.sh` and `Pipeline2_pymol.sh`.
- **20 Nov 2024**: R analysis split into two variants: `Pipeline1_limma.R` (trend regression) and `Pipiline1_ttest.R` (t-test). The `Pipiline` filename typo was introduced here and is preserved in the legacy directory.
- **25 Nov 2024**: BLASTx replaced by BLASTp plus `seqkit translate`; the main FASTQ script renamed to `pipeline1.sh`.

November had **92 commits**. This month established the protein-aware filtering logic and separated the two statistical methods that v2 now merges into `R/differential_genome.R`.

### 2024-12 — Final preprocessing form and pause

- **3 Dec 2024**: `cutadapt` replaced the `seqkit` regex, giving the final adapter-filtering approach.
- **6 Dec 2024**: last update before the repository paused.

December had **5 commits**.

### 2026-09 — v2 refactor

The `refactor/v2` branch was created to address the engineering risks catalogued in `PROJECT_REVIEW.md`: hard-coded paths, no configuration layer, fragile shell constructs, inconsistent preprocessing between FASTQ and POD5 modes, and missing tests. The refactor preserves the legacy pipelines verbatim under `legacy/` and re-implements the workflow as modular, config-driven scripts with a `bin/minbio` dispatcher.

## Commit counts by month

| Month | Commits | Notes |
|---|---|---|
| 2024-09 | 232 | Initial workflow, HPC adaptation, Pipeline1.5, Analysis.R |
| 2024-10 | 196 | seqkit regex, Porechop removal |
| 2024-11 | 92 | BLASTx→BLASTp, Pipeline2 consolidation, limma/t-test split |
| 2024-12 | 5 | cutadapt finalisation |

The drop from September to December reflects a project that reached a working local form and then paused until the v2 refactor.

## Branches

- `main`: original development history; left untouched by the refactor.
- `refactor/v2`: modular, config-driven rewrite with documentation and tests.
- `origin/Reynold4k-Pipeline1`: a stale remote branch from earlier development. It is noted here for provenance; cleanup is optional.

## Interpretation

The commit history shows an iterative, experiment-driven development pattern:

1. A first end-to-end prototype was built quickly.
2. Preprocessing methods were swapped as the author evaluated Porechop, seqkit regex, and finally cutadapt.
3. The protein-level filter evolved from BLASTx to BLASTp plus translation.
4. A second reference mode (Pipeline1.5) was added when full genome annotation was unavailable.
5. Statistical methods were split into trend and t-test variants.
6. After the workflow stabilised, the v2 refactor addressed reproducibility and maintainability without changing the preserved legacy code.

This history is why the v2 design keeps the legacy scripts intact under `legacy/` and documents every behavioural change in `docs/legacy_workflows.md`.
