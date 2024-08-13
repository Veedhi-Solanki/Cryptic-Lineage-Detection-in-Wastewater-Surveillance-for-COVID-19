import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

def list_ids_in_fasta(aligned_fasta_file):
    """
    Lists all IDs in the aligned FASTA file.

    Parameters:
    aligned_fasta_file (str): Path to the aligned FASTA file.
    """
    ids = []
    for record in SeqIO.parse(aligned_fasta_file, "fasta"):
        ids.append(record.id)
    return ids

def generate_consensus(aligned_fasta_file, exclude_branch_ids):
    """
    Generates a consensus sequence from all branches in the aligned FASTA file except the excluded ones.

    Parameters:
    aligned_fasta_file (str): Path to the aligned FASTA file.
    exclude_branch_ids (list): List of branch IDs to exclude.

    Returns:
    Seq: The consensus sequence.
    """
    sequences = [record for record in SeqIO.parse(aligned_fasta_file, "fasta") if record.id not in exclude_branch_ids]
    alignment = MultipleSeqAlignment(sequences)
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    return consensus

def save_consensus_to_file(consensus, output_file):
    """
    Saves the consensus sequence to a FASTA file.

    Parameters:
    consensus (Seq): The consensus sequence.
    output_file (str): The path to the output FASTA file.
    """
    consensus_record = SeqRecord(consensus, id="consensus", description="Consensus sequence excluding specified branches")
    SeqIO.write(consensus_record, output_file, "fasta")

def main():
    # List available aligned FASTA files
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
