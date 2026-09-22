# Troubleshooting

Practical issues and how to resolve them, drawn from the legacy READMEs and the v2 audit.

## BLAST+ open-files limit

**Symptom**: BLAST+ fails with an error like:

```
Error memory mapping: /path/to/db.25.phr openedFilesCount=251 threadID=0
BLAST Database error: Cannot memory map /path/to/db.25.phr. Number of files opened: 251
```

**Fix**: increase the per-process open-files limit before running the workflow:

```bash
ulimit -n unlimited
# or, if the system rejects that:
ulimit -n 65536
```

You can add this to your job script or shell profile.

## Reference/annotation chromosome-name mismatch

**Symptom**: `scripts/check_env.sh` fails with:

```
chr-prefix mismatch between reference FASTA (chr) and GTF (none)
```

**Cause**: the FASTA uses `chr1`, `chr2`, ... while the GTF uses `1`, `2`, ... (or vice versa). `featureCounts` and `samtools view` will silently fail or assign nothing.

**Fix**: use a FASTA and GTF from the same source, or convert chromosome names consistently. `scripts/check_env.sh` catches this before any heavy computation runs.

## Empty BLASTp filter output

**Symptom**: the workflow fails with:

```
filtered FASTQ is empty after BLASTp filtering
```

**Things to inspect**:

- The adapter motif in `filtering.adapter` matches your library (default `GATCCGAATTCN`).
- `cutadapt` is actually finding the motif (`--discard-untrimmed` keeps only reads with the adapter).
- The protein database exists and was built with `makeblastdb -dbtype prot`.
- Translated sequences are long enough to pass `filtering.min_protein_length`.
- Read IDs were not altered between translation and `seqtk subseq` (v2 normalises IDs with `awk '{print $1}'`).

## Windows paths in WSL

**Symptom**: scripts fail to find files when you use Windows-style paths.

**Fix**: in WSL, refer to drives through `/mnt`:

- `D:\Data\experiment` → `/mnt/d/Data/experiment`
- Use forward slashes (`/`) instead of backslashes.

## NanoPlot runs out of memory

**Symptom**: `NanoPlot` is killed or hangs on very large FASTQ files.

**Fix**: run on a compute node with more RAM, or downsample the FASTQ before QC. On this cluster, large tasks should be submitted via `sbatch` rather than run on the login node.

## Dorado model download

**Symptom**: `dorado basecaller` fails because the model is missing.

**Fix**: Dorado can download models automatically on first use if you have network access:

```bash
dorado download --model <model-name>
```

If you are offline, download the model manually from the ONT documentation and point Dorado to the local path. Dorado is not included in the conda environment because of ONT licensing restrictions.

## PyMOL headless rendering

**Symptom**: `workflow/structure_mapping.sh` writes the `.pml` script but no PNG is produced.

**Cause**: PyMOL is not installed or is not on `PATH`.

**Fix**: install `pymol-open-source` from the conda environment, or run the saved script manually:

```bash
pymol -c -q <gene>_highlight.pml
```

If you see OpenGL/headless errors, use an off-screen Mesa build or run PyMOL inside an Xvfb session.

## General debugging tips

- Run `./bin/minbio check --config <yaml>` first. It validates tools, R packages, and reference files.
- Read the timestamped log messages. Every v2 script prints which files it is creating or skipping.
- Use `bash -x workflow/<script>.sh --config <yaml>` for a full trace.
- Ensure all stages are run in order: `check → preprocess → quantify → analyse → candidates → structure`.
