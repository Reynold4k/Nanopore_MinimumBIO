# PROJECT_REVIEW — Nanopore_MinimumBIO

**Reviewer scope**: full repository archaeology (15 files, ~3,900 lines, 525 commits,
2024-09-17 → 2024-12-06), reconstruction of the intended workflow, technical audit,
and refactoring plan. Compiled on the `refactor/v2` branch before any code change.

---

## 1. What the repository currently does

An end-to-end analysis workflow for **Nanopore-sequenced biopanning experiments**
(T7 phage-display peptide libraries screened against small-molecule baits; code and
commit history reference FKBP1B, VHL, CRBN bait–glue systems, JQ1/FK506/Biotin, and a
93-nt personalised T7-Pep reference).

- **Input**: Oxford Nanopore POD5 (raw signal) or FASTQ, organised into panning-round
  directories `R0/R1/.../Rn`, in paired experimental / control folder trees.
- **Processing chain**: basecalling (Dorado) → constant-region / adapter filtering
  (cutadapt `GATCCGAATTCN`) → NanoPlot QC → translation + BLASTp protein-level
  in-frame filtering → alignment (BWA to hg38 **or** minimap2 to a personalised
  target reference) → quantification (featureCounts gene counts **or**
  `samtools depth` coverage) → differential analysis in R (CPM, growth rate,
  trend/t-test, volcano/PCA/line plots) → candidate-region extraction
  (BAM + GTF → BED → FASTA) → Biopython translation → local alignment against an
  AlphaFold PDB sequence → automated PyMOL highlighting.
- **Output**: per-round QC reports, sorted/indexed BAMs, gene-count matrices,
  differential-analysis figures, candidate hit sequences, translated peptides, and
  PyMOL structure highlights (`.pse` + 300 dpi PNG).

The experimental question it answers: *after multiple rounds of biopanning, which
peptides/genes are enriched relative to control, and where do the enriched candidate
sequences map on the 3-D structure of the target protein?*

## 2. Pipeline1 architecture (genome-reference mode)

| File | Role |
|---|---|
| `Pipeline1/pipeline1.sh` | FASTQ entry: merge → **cutadapt** motif filter (finalised 2024-12-03, `6fac716`; earlier iterations used Porechop, then a seqkit regex) → NanoPlot → `seqtk seq -A \| seqkit translate \| blastp` in-frame filter → `bwa mem` → samtools sort/markdup/index → featureCounts |
| `Pipeline1/pod5_preprocessing_with_genecounts.sh` | POD5 entry: npmerge → dorado basecaller → seqkit replace → Porechop → NanoPlot → BWA → sort/markdup → flagstat → featureCounts (**no BLASTp filter** — inconsistent with `pipeline1.sh`) |
| `Pipeline1/Pipeline1_limma.R` | Per-round count matrices → edgeR CPM → exp−ctrl difference matrix → PCA → top-variance trajectory plot → polynomial trend regression `lm(~ poly(time,1))` slope p-values → growth-rate volcano. **Despite the name and `library(limma)`, no limma modelling (`lmFit`/`eBayes`) is used** (rename artefact, `bceeafa`). |
| `Pipeline1/Pipiline1_ttest.R` | Same framework, per-gene unpaired t-test across rounds. Filename typo (`Pipiline`) preserved from history (`7b6e11e`). |

## 3. Pipeline1.5 architecture (personalised / targeted-reference mode)

For experiments where **no genome-wide annotation exists** — only a small custom
reference (e.g. `T7-Pep_Ref_93nt.fasta`). Gene counting is replaced by **coverage**:

- `pipeline1.5.sh`: merge → seqkit motif regex → *deletes original FASTQs*
  (destructive) → NanoPlot → translate + BLASTp filter → `minimap2 -ax map-ont` →
  `samtools depth` → awk aggregation of per-position coverage → top-1000 loci →
  one BED per hit.
- `Compare.sh`: pairs experimental/control rounds, emits `differential_coverage.txt`.
- `Analysis.R` / `Analysis_withoutidmapping.R`: differential coverage → CPM →
  across-round trajectory + log2FC/log10CPM volcano; the former adds a UniProt
  id-mapping step for gene names.

## 4. Pipeline2 architecture (candidate → structural interpretation)

- `Pipeline2_1.sh` (genome mode): pull exon coordinates for a gene from the GTF →
  `samtools view` to extract overlapping reads → `bedtools bamtobed` →
  `bedtools getfasta` DNA extraction → Biopython translation.
- `Pipeline2_1.5.sh` (targeted mode): per-line Python extraction from the custom
  reference using the BEDs produced by Pipeline1.5 → translation.
- `Pipeline2_pymol.sh`: extract chain sequences from an AlphaFold PDB
  (Biopython `PDBParser`) → `PairwiseAligner` local alignment of the translated
  peptide vs PDB sequence → generated PyMOL script highlights the matched segment →
  `.pse` session + ray-traced PNG.

## 5. How the three pipelines relate

Pipeline1 and Pipeline1.5 are **parallel implementations of the same biological
question under two reference conditions** — genome-wide (featureCounts
quantification) vs personalised/targeted reference (coverage quantification): a
two-reference strategy. Pipeline2 is the shared downstream sequence-to-structure
layer for both. Git history confirms iterative, experiment-driven evolution:

```
2024-09-17  initial upload: porechop_preprocessing, dorado, Pipeline1.R, Pipeline2 a/b
2024-09-18  HPC (katana) module-load versions
2024-09-19  Pipeline2 gains BLAST + AlphaFold/PyMOL concepts
2024-09-22  pod5_preprocessing_with_genecounts
2024-09-27  Pipeline1.5 born (targeted library, coverage approach) + Compare
2024-09-30  Analysis.R, Pipeline2_1.5 (personalised branch of Pipeline2)
2024-10-09  seqkit motif regex introduced (GATCCGAATTC[ACGTN])
2024-10-30  Porechop dropped from main FASTQ path
2024-11-04  BLASTx in-frame check added
2024-11-07/08  .pbs → .sh migration (portability off the HPC)
2024-11-15  Pipeline2 consolidated to Pipeline2_1.sh + Pipeline2_pymol.sh
2024-11-20  R analysis split: limma (trend) vs t-test variants
2024-11-25  BLASTx → BLASTp + seqkit translate; renamed to pipeline1.sh
2024-12-03  cutadapt replaces the seqkit regex (final form)
2024-12-06  last update
```

## 6. Complete end-to-end workflow

```
                 ┌─ POD5 ─→ dorado basecaller ─┐
 Raw Nanopore ───┤                             ├─→ FASTQ per round (R0..Rn)
                 └─ FASTQ (direct input) ──────┘
                           │
     motif/adapter filter (cutadapt) + NanoPlot QC
                           │
     protein-level filter: translate (frame 1) → BLASTp vs human proteome
                           │
        ┌─ personalised reference? ─┐
       yes (1.5)                    no (1)
        ▼                            ▼
  minimap2 → samtools depth    bwa mem → markdup → featureCounts
  → top-1000 coverage loci     → gene count matrix
        │                            │
  Compare.sh (diff coverage)   R: CPM / growth rate / trend|t-test
        │                            │
        └──────► volcano / trajectory / PCA → candidate ranking ◄──┘
                           │
                    Pipeline2 (both modes)
     candidate gene/region → BAM subset → BED → DNA → translation
                           │
     AlphaFold PDB sequence ← local alignment → PyMOL highlight figure
```

## 7. Strongest technical contributions (all traceable to code)

1. **End-to-end workflow design** — raw signal → quantification → statistics →
   candidates → protein structure in one repository; few biopanning analyses extend
   to the structural layer.
2. **Protein-aware read filtering** — translate + BLASTp before alignment to remove
   out-of-frame / reverse-orientation reads; deliberately iterated blastx → blastp
   (`1823d63` → `fc52f8d`).
3. **Two-reference strategy** — genome-reference (featureCounts) and
   targeted-reference (coverage) modes for experiments with or without full
   annotation.
4. **Sequence-to-structure integration (Pipeline2)** — BAM→BED→DNA→protein→PDB local
   alignment→PyMOL rendering, fully scripted.
5. **Round-resolved enrichment statistics** — CPM trajectories across panning rounds,
   growth-rate computation, first/last-round background subtraction.

## 8. Weakest engineering aspects

Hard-coded absolute paths (`/mnt/d/...`, `/srv/scratch/z3546698/...`, `D:/...` mixed);
no configuration layer or CLI; no `set -euo pipefail`; no dependency/version pinning;
no tests; half of the README is a Linux beginner tutorial; `Pipiline` typo; duplicate
shebang (`Pipeline2_pymol.sh`); header comments that no longer match the code
(Porechop mentioned after cutadapt replaced it); fragile `while read` subshells.

## 9. Scientific / technical risk audit

Severity ordered. "Confidence" reflects code-reading certainty; items marked
*verify* need an experimental check by the author.

| # | Location | Current implementation | Potential problem | Why it matters | Recommended correction | Confidence |
|---|----------|------------------------|-------------------|----------------|------------------------|------------|
| R1 | `pipeline1.sh` featureCounts step | counts every `*.bam` in step2 | step2 contains both `*_sorted.bam` and `*_marked.bam` — **both counted**; markdup marks but does not remove duplicates | count matrix has two columns per sample; duplicates retained | count only the marked BAM, add `--ignoreDup` (or `markdup -r`) | high |
| R2 | `pod5_preprocessing_with_genecounts.sh` | `dorado basecaller --device cuda "$merged" "$out_dir"` | dorado syntax is `dorado basecaller <model> <data> > out.bam`; model argument missing, no output redirection, output is BAM not a FASTQ tree | POD5 entry **cannot run as written** | `dorado basecaller "$MODEL" "$POD5" > calls.bam` then `samtools fastq` | high |
| R3 | `pipeline1.sh` + README | UCSC `hg38.fa` (chr-prefixed) + Ensembl GTF (no chr prefix) | **chromosome naming mismatch** → featureCounts may assign nothing; `Pipeline2_1.sh` manually adding `"chr"` corroborates this | directly affects quantification correctness | use matched FASTA+GTF sources; validate seqname style at startup | high |
| R4 | `pipeline1.sh`, `pipeline1.5.sh` | `seqtk subseq ... > all_filtered_sequences.fastq.gz` | uncompressed text written with a `.gz` name (downstream survives only via zlib transparency) | file integrity / portability | pipe through `gzip` or drop the `.gz` suffix | high |
| R5 | `pipeline1.sh` translate filter | `awk 'NR%2==1{h=$0} NR%2==0{if(length($0)>=30) print h; print $0}'` | `print $0` sits **outside** the if — length filter is a no-op | the advertised ≥30-aa filter never happens | put both prints inside the if | high |
| R6 | `pipeline1.5.sh` merge | output written inside the searched tree; no exclusion of prior `all_*` files | re-runs fold the previous merged file into itself | duplicated data, irreproducible | exclude `all_*`, collect file list before writing | high |
| R7 | `pipeline1.5.sh` | `find ... ! -name "all*.fastq.gz" -delete` | **deletes the raw sequencing files** | irreversible data loss | never delete inputs; archive instead | high |
| R8 | `Compare.sh` + `Analysis.R` | header written, then file sorted in place; R treats the **last row as column names** | after `sort -k4,4gr` the header lands among zero-difference rows; with negative differences the last row is data | column misassignment or silent skips | keep header out of the sort; R reads with `header=TRUE` | high |
| R9 | `Pipeline1_limma.R` | `exp_base_path <- ".../Control_LuT"`, `control_base_path <- ".../Exp_LuT"` | experimental/control paths **appear swapped** | volcano direction inverted | *author to verify* against lab records; v2 takes paths from config and echoes them | high (needs author confirmation) |
| R10 | `Pipeline1_limma.R` | significance threshold `pvalue < 0.5`; t-test variant uses 0.05; no multiple-testing correction | very liberal + inconsistent | inflated candidate lists | unify threshold, add BH adjustment, document exploratory nature | high |
| R11 | `pipeline1.sh` | `bwa mem` without long-read preset | defaults target Illumina; ONT reads align poorly | mapping rate/accuracy | `-x ont2d` or minimap2 (as in 1.5) | medium-high |
| R12 | `pipeline1.sh` | featureCounts without `-L` | long-read mode required for ONT | counts may be truncated/missing | add `-L` | medium-high |
| R13 | `pipeline1.sh`, pod5 script | `samtools markdup` on ONT reads | designed for short reads; in biopanning, "duplicates" are often genuine enrichment | may flag biological signal | make optional; document rationale | medium |
| R14 | `pipeline1.5.sh` | `seqkit replace` with a regex spanning 4 FASTQ lines | seqkit replace operates per-record on seq/name — a cross-line pattern likely **never matches (silent no-op)** | motif filtering may never have run | use cutadapt (as Pipeline1); *verify* old outputs | medium (verify on target seqkit version) |
| R15 | `pipeline1.sh` BLASTp filter | `seqkit translate` → BLASTp `$1` → `seqtk subseq` | if translate alters read IDs (frame suffix), subseq matching fails → empty filter output | filtering stage may collapse | normalise IDs, assert non-empty output | medium (verify) |
| R16 | `Pipeline2_1.sh` | `awk '$3=="exon" && $0 ~ gene'` substring match over the whole line; whole extracted read regions translated in frame 1 | substring false positives (e.g. `PHF6` vs `PHF6P1`); arbitrary translation frame for non-exonic read spans | wrong coordinates / frame | match the gene_name attribute exactly; document the frame assumption | medium-high |
| R17 | `Pipeline1_limma.R` trend regression | `true_values[1,]` uses only the first column of each round matrix | other per-round BAM columns (incl. R1 duplicates) silently dropped | trend fit on partial data | model all samples explicitly | medium |
| R18 | `Analysis_withoutidmapping.R` | single-round runs **fabricate** a round R(n−1) with CPM = 0 | invented data points in trajectory plots | misleading visualisation | remove; require ≥2 rounds for trajectories | medium-high |
| R19 | R scripts | `library(DESeq2)`, `library(limma)`, `DelayedArray`, `mgcv`, `MASS` loaded but unused | heavy needless dependencies | reproducibility | trim to used packages; pin versions | high |
| R20 | `Analysis_withoutidmapping.R` | identical `mutate(log2foldchange...)` block duplicated | dead code | maintainability | deduplicate | high |
| R21 | pod5 script | `npmerge` for POD5 merging | non-standard tool; official route is `pod5 merge` (or let dorado read a directory) | reproducibility | `pod5 merge` or direct dorado directory input | medium-high |
| R22 | all | no environment definition; only the pod5 script has `module load` (katana-specific) | not reproducible off that HPC | reproducibility | `envs/environment.yml`, R package pin list, `check_env.sh` | high |

These issues are **documented, not silently fixed**: v2 implements the recommended
corrections, and every behavioural difference vs legacy is recorded in
`docs/legacy_workflows.md`.

## 10. Recommended v2 architecture

Modular Bash + R + Python (no Nextflow/Snakemake — disproportionate for this scale),
driven by a single YAML config and a thin dispatcher:

```
README.md
PROJECT_REVIEW.md / REFACTOR_PLAN.md
docs/{workflow,methodology,project_history,legacy_workflows,beginner_setup,troubleshooting}.md
config/config.example.yaml
workflow/{preprocess_fastq,preprocess_pod5,coverage_workflow,candidate_extraction,structure_mapping}.sh
R/{differential_genome.R,differential_coverage.R}
scripts/{lib/common.sh,check_env.sh,translate_dna.py,pdb_align.py}
envs/{environment.yml,r-packages.txt}
examples/{config.fastq.yaml,config.pod5.yaml,config.targeted.yaml,example_tree.md}
tests/smoke_test.sh
bin/minbio                       # subcommand dispatcher
legacy/{Pipeline1,Pipeline1.5,Pipeline2}/ + README.md
```

## 11–13. Preserve / rewrite / deprecate

- **Preserve**: full legacy tree under `legacy/` (filenames as-is, typos included);
  git history (refactor on `refactor/v2`, `main` untouched); the Linux-beginner
  material (moved to `docs/beginner_setup.md`); Pipeline2's sequence-to-structure
  core logic; experimental context examples.
- **Rewrite**: the configuration layer of every script (paths out, `set -euo
  pipefail`, logging, dependency checks, idempotent steps); the R scripts merged and
  parameterised; the main README; heredoc Python extracted to standalone files.
- **Deprecate / flag**: fabricated-round hack (R18); `npmerge` + broken dorado call
  (R2/R21); the `limma` naming (documented, not claimed); raw-FASTQ deletion (R7);
  the stale `origin/Reynold4k-Pipeline1` branch (noted in docs; cleanup optional).

## 14. What the new README highlights (each traceable to code)

- End-to-end Nanopore biopanning workflow: POD5/FASTQ → candidate prioritisation →
  protein structural interpretation
- Protein-aware read filtering (translation + BLASTp) before genomic analysis
- Two-reference strategy: genome-wide (BWA + featureCounts) and targeted
  (minimap2 + coverage) modes
- Round-resolved enrichment statistics (CPM, growth rate, trend regression / t-test)
  with automated visualisation
- Sequence-to-structure integration: candidate regions → DNA → translation →
  AlphaFold PDB alignment → automated PyMOL rendering
- Iterative, experiment-driven development documented in `docs/project_history.md`
