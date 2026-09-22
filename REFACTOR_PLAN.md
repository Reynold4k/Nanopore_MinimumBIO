# REFACTOR_PLAN — Nanopore_MinimumBIO v2

Branch: `refactor/v2`. `main` is never touched. Legacy code is preserved, not deleted.
Scientific behaviour changes vs legacy are listed explicitly in §4 and documented in
`docs/legacy_workflows.md`; nothing is silently altered.

## 1. Objectives

1. Configuration-driven: no hard-coded paths outside `legacy/` and `examples/`.
2. Modular Bash + R + Python; a thin CLI dispatcher (`bin/minbio`).
3. Reproducible: conda env spec, pinned R package list, dependency checker, smoke test.
4. Professional docs: new README, workflow/methodology/history/troubleshooting docs.
5. Preserve legacy: original pipelines live under `legacy/` verbatim.

## 2. Target tree

```
README.md                        # rewritten (see §7)
PROJECT_REVIEW.md                # done
REFACTOR_PLAN.md                 # this file
bin/minbio                       # dispatcher: preprocess|quantify|analyse|candidates|structure|check
config/config.example.yaml
docs/workflow.md                 # end-to-end Mermaid diagram + stage descriptions
docs/methodology.md              # scientific design decisions (two-reference strategy,
                                 #   protein-aware filtering, round-resolved stats, sequence→structure)
docs/project_history.md          # git-history-based evolution (see PROJECT_REVIEW §5)
docs/legacy_workflows.md         # legacy usage + behaviour-change changelog (§4)
docs/beginner_setup.md           # Linux/R/RStudio beginner content migrated from old READMEs
docs/troubleshooting.md
envs/environment.yml             # conda env (see §5)
envs/r-packages.txt              # R package list + install snippet
examples/config.fastq.yaml       # genome mode from FASTQ
examples/config.pod5.yaml        # genome mode from POD5
examples/config.targeted.yaml    # targeted-reference mode
examples/example_tree.md         # expected R0..Rn directory layout
scripts/lib/common.sh            # logging, die, require_tools, yaml config loader
scripts/check_env.sh             # tool + R package presence check, chr-naming sanity check
scripts/translate_dna.py         # Biopython DNA→protein (was heredoc)
scripts/pdb_align.py             # PDB sequence extraction + local alignment (was heredoc)
tests/smoke_test.sh              # bash -n all scripts; toy-config validation; optional toy run
workflow/preprocess_fastq.sh     # genome mode FASTQ entry (was Pipeline1/pipeline1.sh)
workflow/preprocess_pod5.sh      # POD5 entry (fixed dorado call)
workflow/coverage_workflow.sh    # targeted mode (was pipeline1.5.sh + Compare.sh)
workflow/candidate_extraction.sh # BAM/GTF→BED→DNA→protein (merges Pipeline2_1 + 2_1.5)
workflow/structure_mapping.sh    # PDB alignment + PyMOL rendering (was Pipeline2_pymol.sh)
R/differential_genome.R          # merges Pipeline1_limma.R + Pipiline1_ttest.R (--method trend|ttest)
R/differential_coverage.R        # merges Analysis.R + Analysis_withoutidmapping.R (--id-mapping optional)
legacy/Pipeline1/ legacy/Pipeline1.5/ legacy/Pipeline2/   # verbatim copies (git mv)
legacy/README.md                 # what legacy is, how it maps to v2, known issues pointer
```

## 3. Configuration schema (`config/config.example.yaml`)

```yaml
experiment:
  name: "my_experiment"
  exp_folder: "/path/to/experiment"      # contains round dirs R0,R1,...
  control_folder: "/path/to/control"     # same layout
  round_pattern: "R*"                    # glob for round directories

mode: genome                             # genome | targeted

input:
  type: fastq                            # fastq | pod5
  pod5_dir: ""                           # pod5 mode only
  dorado_model: ""                       # pod5 mode only, e.g. dna_r10.4.1_e8.2_400bps_sup@v4.3.0
  dorado_device: "cuda"

reference:
  fasta: "/path/to/reference.fa"
  annotation: "/path/to/annotation.gtf"  # genome mode only
  id_mapping: ""                         # targeted mode, optional UniProt idmapping TSV

filtering:
  adapter: "GATCCGAATTCN"
  cutadapt_error_rate: 0.01
  min_protein_length: 30                 # aa, applied after translation
  protein_db: "/path/to/human_proteome_db"   # makeblastdb protein db
  blast_evalue: 0.001

alignment:
  threads: 16
  aligner: auto                          # auto: bwa (genome) | minimap2 (targeted)
  bwa_preset: "ont2d"                    # applied as -x when aligner=bwa
  mark_duplicates: true
  remove_duplicates: false               # markdup -r when true

analysis:
  method: trend                          # trend | ttest   (genome mode)
  pvalue_threshold: 0.05
  log10cpm_threshold: 3
  top_genes_label: 40
  top_genes_trajectory: 10

structure:
  genes: ["PHF6"]
  uniprot_id: "Q8IWS0"
  pdb_path: "/path/to/AF-Q8IWS0-F1-model_v4.pdb"
  chr_prefix: auto                       # auto | chr | none — seqname style of the BAM
```

Config parsing: `scripts/lib/common.sh` provides `load_config <file>` which calls
`python3` + PyYAML to export `MB_<SECTION>_<KEY>` shell variables. PyYAML is part of
the conda env. Every script also accepts `--config` and errors clearly if unset.

## 4. Behaviour changes vs legacy (must each be documented in docs/legacy_workflows.md)

| Change | Legacy risk ID | New behaviour |
|---|---|---|
| C1 | R1 | featureCounts counts only the marked BAM; `--ignoreDup` when markdup enabled |
| C2 | R2/R21 | dorado called as `dorado basecaller <model> <pod5_dir> \| samtools fastq`; `pod5 merge` optional; no npmerge |
| C3 | R3 | startup check: compare FASTA seqnames vs GTF seqnames; fail loudly on chr-style mismatch |
| C4 | R4 | filtered FASTQ is genuinely gzipped |
| C5 | R5 | ≥30-aa length filter actually applied (both prints inside the awk if) |
| C6 | R6 | merge excludes `all_*` outputs; re-runs are safe (idempotent skip-if-exists) |
| C7 | R7 | raw FASTQs are never deleted |
| C8 | R8 | `differential_coverage.txt` written with header after sorting; R reads `header=TRUE` |
| C9 | R9 | exp/control folders come only from config; echoed at start of every run |
| C10 | R10 | unified `pvalue_threshold` (default 0.05); BH-adjusted p-value column added |
| C11 | R11/R12 | bwa `-x ont2d` (configurable); featureCounts `-L` |
| C12 | R13 | markdup optional via `alignment.mark_duplicates` |
| C13 | R14 | cutadapt used for motif filtering in targeted mode too (not the seqkit cross-line regex) |
| C14 | R15 | after BLASTp filtering, assert non-empty output; IDs normalised with `awk '{print $1}'` |
| C15 | R16 | GTF gene match against the `gene_name` attribute exactly, not whole-line substring |
| C16 | R17 | trend/t-test use all sample columns per round (mean over columns), documented |
| C17 | R18 | no fabricated rounds; trajectory plot requires ≥2 rounds, else skipped with a message |
| C18 | R19/R20 | R scripts load only used packages: edgeR, ggplot2, ggrepel, dplyr, tidyr, rtracklayer, readr |
| C19 | — | `set -euo pipefail`, `require_tools`, timestamped log file per run |

Preserved as-is (documented, not changed): the motif `GATCCGAATTCN`; BLASTp e-value
1e-3 default; top-1000 coverage loci; growth-rate definition (last−first round,
control-subtracted); the PairwiseAligner scoring in structural alignment; PyMOL
rendering style; the polynomial-trend and t-test methods themselves.

## 5. Environments

`envs/environment.yml` (name: minbio): python=3.11, cutadapt, seqkit, seqtk,
blast, bwa, minimap2, samtools, subread, nanoplot, bedtools, biopython, pod5,
pyyaml, pymol-open-source (note: dorado is ONT-licensed — documented manual install).

`envs/r-packages.txt`: CRAN ggplot2/dplyr/tidyr/readr/ggrepel; Bioconductor
edgeR/rtracklayer; with an `renv`-optional note.

## 6. Commit sequence (each commit one logical step)

1. `docs: add project review and refactor plan` ✅ (PROJECT_REVIEW.md, REFACTOR_PLAN.md)
2. `refactor: preserve legacy pipelines under legacy/`
3. `feat: add configuration layer and shared shell library`
4. `refactor: modularise nanopore preprocessing (fastq + pod5)`
5. `feat: add coverage-based targeted-reference workflow`
6. `refactor: consolidate differential analysis in R`
7. `feat: modularise candidate sequence extraction`
8. `feat: add structural mapping workflow`
9. `feat: add reproducible environment and dependency checks`
10. `feat: add minbio CLI dispatcher`
11. `docs: rewrite project README and add user documentation`
12. `docs: document legacy pipeline evolution and migration notes`
13. `test: add smoke tests for core workflow`

## 7. New README outline

`# Nanopore MinimumBIO` + one-sentence positioning; Overview; Key capabilities
(5–8 bullets, all traceable); Workflow (Mermaid, per PROJECT_REVIEW §6);
Workflow modes (genome vs targeted); Installation (conda env + dorado note);
Quick start (`minbio check` → edit config → `minbio preprocess/quantify/analyse/
candidates/structure`); Input (R0..Rn layout, exp/control); Output; Configuration;
Differential analysis; Candidate discovery; Structural interpretation;
Repository structure; Reproducibility; Limitations (honest: exploratory statistics,
translation-frame assumption, chr-naming caveat, dorado licence);
Citation; Author: Chen Zhu. Beginner Linux content moves to docs/beginner_setup.md.

## 8. Acceptance criteria

- `bash -n` passes for every shell script; `Rscript -e 'parse(...)'` passes for R.
- No absolute `/mnt/`, `/srv/`, `D:/` paths outside `legacy/` (grep check).
- `tests/smoke_test.sh` runs and passes without needing the heavy tools installed
  (syntax + config validation level), and detects+reports missing tools clearly.
- `bin/minbio --help` lists subcommands; each maps to the right workflow script.
- Every claim in README's "Key capabilities" points at a real script.
- docs/legacy_workflows.md contains the full §4 changelog table.
