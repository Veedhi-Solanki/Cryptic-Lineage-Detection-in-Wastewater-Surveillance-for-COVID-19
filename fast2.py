import os
import subprocess
from Bio import Phylo, SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def filter_short_sequences(fasta_file, min_length=151):
    filtered_fasta_file = f"{fasta_file}_filtered.fasta"
    with open(fasta_file, "r") as input_handle, open(filtered_fasta_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record.seq) >= min_length:
                SeqIO.write(record, output_handle, "fasta")
    return filtered_fasta_file

def cluster_sequences(fasta_file, threshold=0.99):
    clustered_fasta_file = f"{fasta_file}_clustered.fasta"
    print(f"Clustering sequences in {fasta_file} with a threshold of {threshold}...")
    subprocess.run(f"cd-hit -i {fasta_file} -o {clustered_fasta_file} -c {threshold} -n 5 -d 0 -T 8", shell=True, check=True)
    return clustered_fasta_file

def build_tree(sample_id, bam_type='sorted'):
    bam_file = f"{sample_id}_{bam_type}.bam"
    fasta_file = f"{sample_id}.fasta"
    phylogenetic_tree_file = f"{fasta_file}.treefile"

    if not os.path.isfile(bam_file):
        print(f"BAM file {bam_file} does not exist.")
        return

    print(f"Converting {bam_file} to FASTA...")
    subprocess.run(f"samtools fasta {bam_file} > {fasta_file}", shell=True, check=True)

    print(f"Filtering short sequences from {fasta_file}...")
    filtered_fasta_file = filter_short_sequences(fasta_file)

    print(f"Clustering sequences to reduce redundancy...")
    clustered_fasta_file = cluster_sequences(filtered_fasta_file, threshold=0.90)  # Adjust threshold here

    print(f"Building phylogenetic tree for {clustered_fasta_file} with FastTree...")
    subprocess.run(f"FastTree -nt {clustered_fasta_file} > {phylogenetic_tree_file}", shell=True, check=True)

    return phylogenetic_tree_file

def save_tree_as_image(tree_file_path):
    print(f"Reading and saving the phylogenetic tree from {tree_file_path} as an image...")

    # Read the tree file
    tree = Phylo.read(tree_file_path, "newick")
    
    # Create a plot with a larger figure size
    fig = plt.figure(figsize=(20, 20))  # Adjusted size for better readability
    ax = fig.add_subplot(1, 1, 1)
    
    # Draw the tree with adjustments to labels
    Phylo.draw(tree, do_show=False, axes=ax)
    
    # Rotate labels and adjust spacing
    plt.xticks(rotation=90, ha='center')
    plt.tight_layout()
    
    # Save the plot as a PNG file
    image_file_path = f"{tree_file_path}_phylogenetic_tree.png"
    plt.savefig(image_file_path, bbox_inches='tight')
    plt.close()
    
    print(f"Phylogenetic tree saved as {image_file_path}")

if __name__ == "__main__":
    sample_id = input("Enter the sample ID: ")
    tree_file = build_tree(sample_id)
    if tree_file:
        save_tree_as_image(tree_file)
