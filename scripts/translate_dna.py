#!/usr/bin/env python3
"""Translate DNA FASTA records into protein sequences.

Standalone replacement for the legacy Biopython heredoc used in Pipeline2.
Pads sequences to a multiple of three with 'N' and translates with
``to_stop=False`` so that internal stop codons become 'X'.
"""

import argparse
import sys

from Bio.Seq import Seq
from Bio import SeqIO


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Translate DNA FASTA records into protein FASTA."
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input DNA FASTA file."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output protein FASTA file."
    )
    args = parser.parse_args()

    written = 0
    with open(args.output, "w") as f_out:
        for record in SeqIO.parse(args.input, "fasta"):
            dna_seq = Seq(str(record.seq))
            print(
                f"[translate_dna] {record.id}: length {len(dna_seq)}",
                file=sys.stderr,
                flush=True,
            )

            remainder = len(dna_seq) % 3
            if remainder != 0:
                pad_length = 3 - remainder
                dna_seq = dna_seq + "N" * pad_length
                print(
                    f"[translate_dna] {record.id}: padded to {len(dna_seq)}",
                    file=sys.stderr,
                    flush=True,
                )

            protein_seq = dna_seq.translate(to_stop=False)
            f_out.write(f">{record.id}\n{protein_seq}\n")
            written += 1

    print(
        f"[translate_dna] Wrote {written} protein records to {args.output}",
        file=sys.stderr,
        flush=True,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
