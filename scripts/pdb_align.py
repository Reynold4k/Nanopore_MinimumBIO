#!/usr/bin/env python3
"""Extract PDB chain sequences and local-align a query protein.

Standalone replacement for the Python heredoc blocks in
legacy/Pipeline2/Pipeline2_pymol.sh. Uses Biopython's PDBParser and
PairwiseAligner with the legacy scoring parameters.
"""

import argparse
import sys

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBParser

AA3TO1 = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}


def extract_chain_sequences(pdb_path: str):
    """Yield (chain_id, sequence) tuples from the first model."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("Protein", pdb_path)
    model = structure[0]
    for chain in model:
        residues = [r for r in chain if r.get_resname() in AA3TO1]
        seq = "".join(AA3TO1[r.get_resname()] for r in residues)
        yield chain.id, seq


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Local-align a query protein against PDB chain sequences."
    )
    parser.add_argument("--pdb", required=True, help="Input PDB file.")
    parser.add_argument(
        "--query", required=True, help="Query protein FASTA file."
    )
    args = parser.parse_args()

    records = list(SeqIO.parse(args.query, "fasta"))
    if not records:
        print("No sequences in query FASTA", file=sys.stderr)
        return 1

    protein_seq = records[0].seq
    print(f"[pdb_align] Query: {records[0].id} ({len(protein_seq)} aa)", file=sys.stderr)

    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    best_score = 0
    best_chain = None
    best_position = (0, 0)

    for chain_id, pdb_seq in extract_chain_sequences(args.pdb):
        print(
            f"[pdb_align] Chain {chain_id}: {len(pdb_seq)} aa",
            file=sys.stderr,
            flush=True,
        )
        for alignment in aligner.align(protein_seq, pdb_seq):
            aligned_len = alignment.aligned[0][0][1] - alignment.aligned[0][0][0]
            if alignment.score > best_score and aligned_len > 5:
                best_score = alignment.score
                best_chain = chain_id
                best_position = alignment.aligned[1][0]

    if best_chain is None:
        print("No alignment found", file=sys.stderr)
        return 1

    start_index, end_index = best_position
    end_index -= 1  # Convert to zero-based inclusive.
    print(f"{best_chain} {start_index + 1} {end_index + 1}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
