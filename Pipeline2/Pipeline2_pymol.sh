#!/bin/bash

# This script facilitates the analysis and visualization of genomic and protein data by integrating various bioinformatics tools
# and techniques. It extracts specific genomic regions, translates DNA sequences into proteins, matches these translated proteins to 
# structural protein data obtained from PDB files, and generates visualizations using PyMOL. By aligning protein sequences obtained 
# from high-coverage DNA sequences to reference structures in PDB format, the script identifies the best matching regions,
# highlighting them for further analysis. The workflow includes sequence extraction, translation, local sequence alignment,
# and structural visualization, ultimately aiding in visual insights into how specific genomic variants or regions might affect 
# protein structure and function. This automated approach helps streamline molecular biology and genomics research by simplifying 
# sequence-to-structure comparisons and generating publication-ready images.


# What you need to modify:

# 1.Change the line 24 VISUAL_DIR to your/path/to/exp/potential_hit/visualization
# 2.Change the line 27 GENOME_FASTA to your/path/to/reference
# 3.Change the line 28 BED_FILE to your/path/to/exp/potential_hit/gene_Hit_all_trimmed_sorted_merged.bed file
# 4.Change the line 32 UNIPROT_ID to your gene's uniprot id
# 5.Do not change the last part of the sample directory! Change the line 35 PDB_PATH to your/path/to/alphafold database/UP000005640_9606_HUMAN_v4/AF-${UNIPROT_ID}-F1-model_v4.pdb



#!/bin/bash

# Define the output directory
VISUAL_DIR="/mnt/d/A_FKBP1B/WDB001/YB/241004/potential_hit/visualization"  # Directory for visualizations
mkdir -p "$VISUAL_DIR"  # Ensure the directory exists

# Define key file paths and directories
GENOME_FASTA="/mnt/d/hg38/hg38.fa" # Reference FASTA file
BED_FILE="/mnt/d/A_FKBP1B/WDB001/YB/241004/potential_hit/PHF6_Hit_all_filtered_sequences_sorted_sorted.bed"  # BED file

# PDB and sequence paths
UNIPROT_ID="Q8IWS0"
PDB_PATH="/mnt/d/hg38/alphafold/AF-${UNIPROT_ID}-F1-model_v4.pdb"

FILENAME=$(basename "$BED_FILE")
GENE_NAME=${FILENAME%%_*}
BED_DIR=$(dirname "$BED_FILE")

dna_sequences_output="$BED_DIR/high_coverage_sequences.fa"
translated_proteins_output="$BED_DIR/translated_proteins.fa"
pdb_fasta_file="$BED_DIR/pdb_sequence_${UNIPROT_ID}.fasta"

# Extract protein sequence from translated_proteins.fa
protein_sequence=$(awk '/^>/ {if (seqlen) exit} {if (!/^>/) printf $0}' "$translated_proteins_output")
echo "Extracted Protein Sequence (before cleanup): '$protein_sequence'"

# Clean protein sequence, allowing 'X'
protein_sequence=$(echo "$protein_sequence" | tr -d '[:space:]' | tr -d '*')
echo "Cleaned Protein Sequence: '$protein_sequence'"

# Skip validation and proceed to prepare the PDB sequence in FASTA format using Biopython
python3 -c "
from Bio import SeqIO
from Bio.PDB import PDBParser

# Load PDB and extract sequences
pdb_file = '$PDB_PATH'
pdb_fasta = '$pdb_fasta_file'

parser = PDBParser(QUIET=True)
structure = parser.get_structure('Protein', pdb_file)
model = structure[0]

# Mapping three-letter codes to one-letter
aa3to1 = {'CYS':'C', 'ASP':'D', 'SER':'S', 'GLN':'Q', 'LYS':'K', 
          'ILE':'I', 'PRO':'P', 'THR':'T', 'PHE':'F', 'ASN':'N', 
          'GLY':'G', 'HIS':'H', 'LEU':'L', 'ARG':'R', 'TRP':'W', 
          'ALA':'A', 'VAL':'V', 'GLU':'E', 'TYR':'Y', 'MET':'M'}

sequences = []

for chain in model:
    residues = [residue for residue in chain if residue.get_resname() in aa3to1]
    seq = ''.join(aa3to1[residue.get_resname()] for residue in residues)
    sequences.append((f'Chain_{chain.id}', seq))

with open(pdb_fasta, 'w') as f:
    for chain_id, sequence in sequences:
        f.write(f'>{chain_id}\n{sequence}\n')
"

# Verify that pdb_fasta_file is created
if [[ ! -f $pdb_fasta_file ]]; then
    echo "PDB FASTA file creation failed."
    exit 1
fi

alignment_results=$(python3 - <<END
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

# Load the sequences
protein_seq = Seq("$protein_sequence")

with open('$pdb_fasta_file', 'r') as f:
    pdb_records = list(SeqIO.parse(f, 'fasta'))

# Initialize aligner
aligner = PairwiseAligner()
aligner.mode = 'local'  # Use local alignment

# Adjust alignment parameters for longer matches
aligner.open_gap_score = -2  # Adjust gap open penalty
aligner.extend_gap_score = -0.5  # Adjust gap extension penalty
aligner.match_score = 1  # Set higher match score
aligner.mismatch_score = 0  # Set mismatch penalty

best_score = 0
best_alignment = None
best_chain = None
best_position = (0, 0)

# Evaluate each chain in the PDB
for record in pdb_records:
    pdb_seq = record.seq
    alignments = aligner.align(protein_seq, pdb_seq)

    # Check all alignments to find the best one
    for alignment in alignments:
        if alignment.score > best_score and (alignment.aligned[0][0][1] - alignment.aligned[0][0][0]) > 5: # Example threshold for length
            best_score = alignment.score
            best_alignment = alignment
            best_chain = record.id
            best_position = alignment.aligned[1][0]  # This is (start, end)

start_index, end_index = best_position
end_index -= 1  # Convert to zero-based inclusive

print(f'{best_chain} {start_index + 1} {end_index + 1}')
END
)

# Capture the best alignment details
read chain_id start_index end_index <<< "$alignment_results"

echo "Chain ID: $chain_id"
echo "Start Index: $start_index"
echo "End Index: $end_index"

chain_id="${chain_id#Chain_}"

if [[ -n $chain_id && -n $start_index && -n $end_index ]]; then
    echo "Creating PyMOL script with chain $chain_id and resi $start_index-$end_index"
    SCRIPT=$(cat <<-EOF
    load $PDB_PATH, structure
    show cartoon, structure
    select highlight_region, chain $chain_id and resi $start_index-$end_index
    color gray, structure
    color red, highlight_region
    save $VISUAL_DIR/highlighted_structure.pse, format=pse
    png $VISUAL_DIR/${GENE_NAME}_highlighted_structure.png, width=800, height=600, dpi=300, ray=1
    quit
EOF
    )
    pymol -c -d "${SCRIPT}"
else
    echo "Variable issue: ensure chain_id and indices are correctly defined."
fi
