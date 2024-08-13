import os
import subprocess
from Bio import Phylo, SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def filter_short_sequences(fasta_file, min_length=151):
    # Function to filter out short sequences from a FASTA file
    filtered_fasta_file = f"{fasta_file}_filtered.fasta"
    with open(fasta_file, "r") as input_handle, open(filtered_fasta_file, "w") as output_handle:
        # Parse the FASTA file and only write sequences longer than the specified minimum length
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record.seq) >= min_length:
                SeqIO.write(record, output_handle, "fasta")
    return filtered_fasta_file

def cluster_sequences(fasta_file, threshold=0.99):
    # Function to cluster sequences using CD-HIT to reduce redundancy
    clustered_fasta_file = f"{fasta_file}_clustered.fasta"
    print(f"Clustering sequences in {fasta_file} with a threshold of {threshold}...")
    subprocess.run(f"cd-hit -i {fasta_file} -o {clustered_fasta_file} -c {threshold} -n 5 -d 0 -T 8", shell=True, check=True)
    return clustered_fasta_file

def build_tree(sample_id, bam_type='sorted'):
    # Main function to build a phylogenetic tree from a BAM file
    bam_file = f"{sample_id}_{bam_type}.bam"
    fasta_file = f"{sample_id}.fasta"
    phylogenetic_tree_file = f"{fasta_file}.treefile"

    # Check if the BAM file exists
    if not os.path.isfile(bam_file):
        print(f"BAM file {bam_file} does not exist.")
        return

    # Convert BAM file to FASTA format using samtools
    print(f"Converting {bam_file} to FASTA...")
    subprocess.run(f"samtools fasta {bam_file} > {fasta_file}", shell=True, check=True)

    # Filter out short sequences from the generated FASTA file
    print(f"Filtering short sequences from {fasta_file}...")
    filtered_fasta_file = filter_short_sequences(fasta_file)

    # Cluster the sequences to reduce redundancy using CD-HIT
    print(f"Clustering sequences to reduce redundancy...")
    clustered_fasta_file = cluster_sequences(filtered_fasta_file, threshold=0.90)  # Adjust threshold here

    # Build a phylogenetic tree using FastTree with the clustered FASTA file
    print(f"Building phylogenetic tree for {clustered_fasta_file} with FastTree...")
    subprocess.run(f"FastTree -nt {clustered_fasta_file} > {phylogenetic_tree_file}", shell=True, check=True)

    return phylogenetic_tree_file

def save_tree_as_image(tree_file_path):
    # Function to save the phylogenetic tree as an image
    print(f"Reading and saving the phylogenetic tree from {tree_file_path} as an image...")

    # Read the tree file in Newick format
    tree = Phylo.read(tree_file_path, "newick")
    
    # Create a plot with a larger figure size for better readability
    fig = plt.figure(figsize=(20, 20))  # Adjusted size for better readability
    ax = fig.add_subplot(1, 1, 1)
    
    # Draw the phylogenetic tree without displaying it (do_show=False)
    Phylo.draw(tree, do_show=False, axes=ax)
    
    # Rotate labels and adjust spacing for readability
    plt.xticks(rotation=90, ha='center')
    plt.tight_layout()
    
    # Save the plot as a PNG image with tight bounding box
    image_file_path = f"{tree_file_path}_phylogenetic_tree.png"
    plt.savefig(image_file_path, bbox_inches='tight')
    plt.close()
    
    print(f"Phylogenetic tree saved as {image_file_path}")

if __name__ == "__main__":
    # Prompt the user for the sample ID
    sample_id = input("Enter the sample ID: ")
    # Build the phylogenetic tree for the provided sample ID
    tree_file = build_tree(sample_id)
    # If the tree was successfully built, save it as an image
    if tree_file:
        save_tree_as_image(tree_file)
