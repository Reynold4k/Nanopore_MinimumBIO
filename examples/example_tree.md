# Expected input directory layout

The experiment and control folders follow the same round-based structure.
Each round directory is named `R<n>` (e.g. `R0`, `R1`, `R2`) and contains the
FASTQ files produced for that selection round.

## Genome mode from FASTQ (`examples/config.fastq.yaml`)

```
/path/to/experiment/
├── R0/
│   ├── round0_sample1.fastq.gz
│   └── round0_sample2.fastq.gz
├── R1/
│   └── round1.fastq.gz
├── R2/
│   └── round2.fastq.gz
└── ...

/path/to/control/
├── R0/
│   └── control_round0.fastq.gz
├── R1/
│   └── control_round1.fastq.gz
└── ...
```

The `round_pattern` config key (default `R*`) is used to discover round
directories. Round names must match between the experiment and control folders
so that pairs can be compared.

## Genome mode from POD5 (`examples/config.pod5.yaml`)

```
/path/to/pod5/
├── batch0.pod5
├── batch1.pod5
└── ...

/path/to/experiment/
├── R0/
├── R1/
└── ...
```

Dorado basecalls the POD5 files and writes FASTQ into the round directories.
After basecalling, the layout matches the FASTQ genome mode above.

## Targeted mode (`examples/config.targeted.yaml`)

Targeted mode uses the same round-based FASTQ layout as genome FASTQ mode, but
aligns against a smaller, targeted reference (e.g. a synthetic library or
amplicon) instead of a full genome. The optional UniProt id-mapping TSV can be
used to annotate coverage results.
