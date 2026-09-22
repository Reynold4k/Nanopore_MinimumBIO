# Methodology

This document explains the scientific design decisions behind Nanopore MinimumBIO and the assumptions that shape the analysis.

## Experimental context

The workflow was designed for **T7 phage-display peptide libraries** screened against small-molecule baits. Each panning round enriches peptides that bind the bait. Sequencing those rounds with Oxford Nanopore gives long reads that still contain the library constant region and the variable peptide insert. The pipeline asks: *which peptides or genes become enriched across rounds relative to a control, and where do the enriched sequences map on the target protein structure?*

Reads are organised as paired experimental and control folder trees, each containing round directories (`R0`, `R1`, ... `Rn`). Round `R0` is typically the unselected input library; later rounds are progressively enriched.

## Protein-aware filtering

Before any alignment, reads are filtered at the protein level:

1. The constant-region motif `GATCCGAATTCN` is removed with `cutadapt`.
2. Surviving reads are translated in reading frame 1 with `seqkit translate`.
3. Translated peptides ≥ `filtering.min_protein_length` amino acids are screened with `blastp` against the configured protein database.
4. Only reads whose peptide has a significant BLAST hit are kept for alignment.

This removes out-of-frame, reverse-orientation, and non-library reads before they can bias genomic quantification. The motif `GATCCGAATTCN` is treated as the library constant region and is preserved from the original pipeline.

## Two-reference strategy

The same biological question is answered under two reference regimes:

- **Genome mode**: align to a full genome (`bwa mem`), count reads per gene with `featureCounts`, and run differential gene-expression statistics. Used when a genome-wide reference and GTF annotation are available.
- **Targeted mode**: align to a small custom reference such as a synthetic library or amplicon (`minimap2`), compute per-position coverage with `samtools depth`, and compare coverage between experiment and control. Used when no full annotation exists.

This dual design lets the same repository handle both standard human-genome screens and custom T7-Pep libraries.

## Round-resolved enrichment statistics

### Count normalisation

Genome-mode gene counts are converted to counts per million (CPM) with edgeR. Targeted-mode coverage values are also normalised to CPM within each round.

### Growth rate and background subtraction

The pipeline computes two related quantities:

- **Per-round difference**: mean experimental CPM minus mean control CPM for the round.
- **Growth rate**: the legacy definition is `(last round difference) - (first round difference)`. Positive values indicate enrichment over the course of panning.

These definitions are preserved from the legacy scripts because they match the lab's existing interpretation.

### Differential methods

Genome mode supports two methods, chosen with `analysis.method`:

- **Trend**: a first-degree polynomial regression of the per-round experiment–control difference against round index. The slope p-value ranks genes whose enrichment changes across rounds.
- **t-test**: an unpaired t-test comparing experimental and control CPM vectors per gene across rounds.

Both methods treat panning rounds as ordered time points rather than as independent biological replicates. P-values are therefore best used for ranking candidates, not for formal statistical inference.

## Sequence-to-structure rationale

After differential ranking, candidate regions are translated and mapped to a 3-D structure:

1. For genome mode, exon coordinates for the candidate gene are pulled from the GTF and reads overlapping those coordinates are extracted from the BAM.
2. For targeted mode, the top-coverage BEDs are used directly.
3. Extracted BED intervals are converted back to DNA with `bedtools getfasta`.
4. DNA is translated in frame 1 with `scripts/translate_dna.py`.
5. Translated peptides are aligned locally against the AlphaFold PDB sequence using a Biopython `PairwiseAligner`.
6. The best alignment is used to generate a PyMOL script that highlights the matched residue range.

This closes the loop from sequencing reads to a structural hypothesis about where an enriched peptide binds.

## Explicit caveats

- **Exploratory statistics**: the workflow does not require biological replicates and treats rounds as repeated measures. P-values are for ranking, not for rigorous hypothesis testing.
- **No replicates assumption**: because each round is a single selection step, there is no within-round replication. Results should be validated experimentally.
- **Frame-1 translation assumption**: protein-aware filtering and candidate translation both use reading frame 1. Peptides enriched in other frames will be missed.
- **Chromosome naming**: the FASTA and GTF must share the same chromosome style. A mismatch is caught by `scripts/check_env.sh`.
- **Targeted mode coverage aggregation**: coverage is summed per reference sequence and the top 1000 intervals are reported. This is a coarse, exploratory view of enrichment.
- **Structural alignment**: local alignment finds the best matching segment but does not model binding affinity or confirm a biological interaction.
