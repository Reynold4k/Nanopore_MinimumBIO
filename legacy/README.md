# Legacy Pipelines

This directory preserves the original Nanopore MinimumBIO pipeline implementations
verbatim for provenance and reproducibility.

- `legacy/Pipeline1/` — original genome-mode preprocessing, gene counting and
differential analysis (`pipeline1.sh`, `Pipeline1_limma.R`, `Pipiline1_ttest.R`,
`pod5_preprocessing_with_genecounts.sh`).
- `legacy/Pipeline1.5/` — original targeted-reference coverage workflow
(`pipeline1.5.sh`, `Compare.sh`, `Analysis.R`, `Analysis_withoutidmapping.R`).
- `legacy/Pipeline2/` — original candidate sequence extraction and structural
mapping (`Pipeline2_1.sh`, `Pipeline2_1.5.sh`, `Pipeline2_pymol.sh`).

These files are kept exactly as they were in the pre-v2 repository, including the
historical filename typo `Pipiline1_ttest.R`.

For the refactored v2 workflows, see:

- `../REFACTOR_PLAN.md` §4 for a table of scientific and behaviour changes between
  the legacy pipelines and the new modular implementation.
- `../docs/legacy_workflows.md` for legacy usage instructions and the full
  changelog.

Do not modify the contents of `legacy/*`; any bug fixes or improvements belong in
`workflow/` or `scripts/`.
