import os
import subprocess
from Bio import Phylo
import matplotlib.pyplot as plt

def build_tree(sample_id, bam_type='sorted'):
    # Define file names based on sample ID and BAM file type
    bam_file = f"{sample_id}_{bam_type}.bam"
    fastq_file = f"{sample_id}.fastq"
    fasta_file = f"{sample_id}.fasta"
    aligned_fasta_file = f"{sample_id}_aligned.fasta"
    phylogenetic_tree_file = f"{aligned_fasta_file}.treefile"

    # Check if the BAM file exists
    if not os.path.isfile(bam_file):
        print(f"BAM file {bam_file} does not exist.")
        return

    # Convert BAM file to FASTA format using samtools and seqtk
    print(f"Converting {bam_file} to FASTA...")
    subprocess.run(f"samtools bam2fq {bam_file} | seqtk seq -A > {fasta_file}", shell=True, check=True)

    # Align the FASTA sequences using MAFFT
    print(f"Aligning sequences of {fasta_file} with MAFFT...")
    subprocess.run(f"mafft --auto {fasta_file} > {aligned_fasta_file}", shell=True, check=True)

    # Build a phylogenetic tree using IQ-TREE with the aligned FASTA file
    print(f"Building phylogenetic tree for {aligned_fasta_file} with IQ-TREE...")
    subprocess.run(f"iqtree -s {aligned_fasta_file} -nt AUTO -m GTR -redo", shell=True, check=True)

    # Read and display the phylogenetic tree
    print(f"Reading and displaying the phylogenetic tree for {sample_id}...")
    tree = Phylo.read(phylogenetic_tree_file, "newick")
    Phylo.draw(tree)

    # Save the displayed phylogenetic tree as a PNG image
    plt.savefig(f"{sample_id}_phylogenetic_tree.png")
    plt.show()

if __name__ == "__main__":
    # Prompt the user for the sample ID and BAM file type
    sample_id = input("Enter the sample ID: ")
    bam_type = input("Enter the BAM file type (aligned, sorted, mapped): ")
    build_tree(sample_id, bam_type)
