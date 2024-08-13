import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def extract_sequence(aligned_fasta_file, branch_id):
    """
    Extracts the sequence of the specified branch from the aligned FASTA file.

    Parameters:
    aligned_fasta_file (str): Path to the aligned FASTA file.
    branch_id (str): Identifier of the branch to extract.

    Returns:
    SeqRecord: The sequence record of the specified branch.
    """
    for record in SeqIO.parse(aligned_fasta_file, "fasta"):
        if record.id == branch_id:
            return record
    return None

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

def save_sequence_to_file(sequence, output_file):
    """
    Saves the sequence record to a FASTA file.

    Parameters:
    sequence (SeqRecord): The sequence record to save.
    output_file (str): The path to the output FASTA file.
    """
    SeqIO.write(sequence, output_file, "fasta")

def main():
    # List available aligned FASTA files
    print("Available aligned FASTA files in the current directory:")
    for file in os.listdir('.'):
        if file.endswith('_aligned.fasta'):
            print(file)

    # Input ID number and node number
    id_number = input("Enter the ID number (e.g., SRR16241857): ")
    node_number = input("Enter the node number (e.g., 92-32): ")

    # Construct the aligned FASTA file name
    aligned_fasta_file = f"{id_number}_aligned.fasta"

    # Use the node number directly as branch ID
    branch_id = node_number

    # Check if the file exists
    if not os.path.isfile(aligned_fasta_file):
        print(f"Aligned FASTA file {aligned_fasta_file} does not exist.")
        return

    # List all IDs in the aligned FASTA file
    ids = list_ids_in_fasta(aligned_fasta_file)
    print("IDs in the aligned FASTA file:")
    print(ids)

    # Extract and print the sequence of the long branch
    long_branch_seq = extract_sequence(aligned_fasta_file, branch_id)
    if long_branch_seq:
        print(f"Sequence of the long branch ({branch_id}):")
        print(long_branch_seq.seq)

        # Save the sequence to a new FASTA file
        output_file = f"{id_number}_long_branch_{branch_id}.fasta"
        save_sequence_to_file(long_branch_seq, output_file)
        print(f"Sequence saved to {output_file}")
    else:
        print(f"Branch ID {branch_id} not found in the aligned FASTA file.")

if __name__ == "__main__":
    main()
