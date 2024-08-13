import os
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio import motifs

def list_ids_in_fasta(fasta_file):
    """List all sequence IDs in a FASTA file."""
    ids = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids.append(record.id)
    return ids

def generate_consensus(fasta_file, exclude_branch_ids):
    """Generate a consensus sequence from an aligned FASTA file, excluding specified branches."""
    alignment = AlignIO.read(fasta_file, "fasta")
    included_records = [record for record in alignment if record.id not in exclude_branch_ids]
    included_alignment = MultipleSeqAlignment(included_records)
    
    motif = motifs.create(included_alignment)
    consensus_seq = motif.consensus
    
    return consensus_seq

def save_consensus_to_file(consensus_seq, output_file):
    """Save the consensus sequence to a FASTA file."""
    with open(output_file, "w") as handle:
        handle.write(f">Consensus\n{consensus_seq}\n")

def main():
    print("Available aligned FASTA files in the current directory:")
    for file in os.listdir('.'):
        if file.endswith('_aligned.fasta'):
            print(file)

    # Input ID number
    id_number = input("Enter the ID number (e.g., SRR16241857): ")

    # Construct the aligned FASTA file name
    aligned_fasta_file = f"{id_number}_aligned.fasta"

    # Check if the file exists
    if not os.path.isfile(aligned_fasta_file):
        print(f"Aligned FASTA file {aligned_fasta_file} does not exist.")
        return

    # List all IDs in the aligned FASTA file
    ids = list_ids_in_fasta(aligned_fasta_file)
    print("IDs in the aligned FASTA file:")
    print(ids)

    # Input IDs to exclude
    exclude_branch_ids = input("Enter branch IDs to exclude separated by commas (e.g., 92-32,1-33232): ").split(',')

    # Generate consensus sequence for all branches except the excluded ones
    consensus_seq = generate_consensus(aligned_fasta_file, exclude_branch_ids)
    print("Consensus sequence excluding specified branches:")
    print(consensus_seq)

    # Save the consensus sequence to a file
    output_file = f"{id_number}_consensus.fasta"
    save_consensus_to_file(consensus_seq, output_file)
    print(f"Consensus sequence saved to {output_file}")

if __name__ == "__main__":
    main()
