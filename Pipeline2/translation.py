from Bio import SeqIO
import os

def translate_contigs_to_protein(fasta_file, protein_file):
    with open(protein_file, "w") as output:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Assume reading frame is chosen to encompass longest ORF or by further intelligence
            for frame in range(3):  # Attempt translations in all frames
                dna_seq = record.seq[frame:]
                protein_seq = dna_seq.translate(to_stop=True)
                if len(protein_seq) > 30:  # Arbitrary threshold for meaningful protein length
                    output.write(f">{record.id}_frame_{frame}\n{protein_seq}\n")
                    break

# Example usage
contig_fasta_file = "/srv/scratch/z3546698/true/Small_Molecule/JQ1/T7MB-1/231104/potential_hit/FASTQ/canu_out/corrected.contigs.fasta"
output_protein_file = "/srv/scratch/z3546698/true/Small_Molecule/JQ1/T7MB-1/231104/potential_hit/FASTQ/canu_out/protein_sequences_long.fasta"
translate_contigs_to_protein(contig_fasta_file, output_protein_file)
