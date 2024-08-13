import os
from Bio import SeqIO

def list_ids_in_fasta(fasta_file):
    """List all sequence IDs in a FASTA file."""
    ids = []
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            ids.append(record.id)
    return ids

def extract_sequence(fasta_file, sequence_id):
    """Extract a sequence from a FASTA file by its ID."""
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == sequence_id:
                return record
    return None

def save_sequence_to_file(sequence, output_file):
    """Save a single sequence to a FASTA file."""
    with open(output_file, "w") as handle:
        SeqIO.write(sequence, handle, "fasta")

def main():
    # Input ID number and sequence ID
    id_number = input("Enter the ID number (e.g., As20221013_S20_L001): ")
    sequence_id = input("Enter the sequence ID to extract (e.g., 02023-26): ")

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

    # Extract and print the sequence with the specified ID
    sequence = extract_sequence(aligned_fasta_file, sequence_id)
    if sequence:
        print(f"Sequence of {sequence_id}:")
        print(sequence.seq)

        # Save the sequence to a new FASTA file
        output_file = f"{id_number}_sequence_{sequence_id}.fasta"
        save_sequence_to_file(sequence, output_file)
        print(f"Sequence saved to {output_file}")
    else:
        print(f"Sequence ID {sequence_id} not found in the aligned FASTA file.")

if __name__ == "__main__":
    main()
