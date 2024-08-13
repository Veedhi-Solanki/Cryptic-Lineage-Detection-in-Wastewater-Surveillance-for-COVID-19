# Cryptic Lineage Detection in Wastewater COVID-19 Surveillance

This repository contains a series of scripts designed for the detection of cryptic SARS-CoV-2 lineages in wastewater samples. The scripts are organized and numbered according to the order in which they should be executed for both Ontario and NYC workflows. 

**Note:** The code in the files is identical for both Ontario and NYC workflows, except for the `build_tree` scripts, where different software is used: FastTree for Ontario and IQ-TREE for NYC.

## Project Structure

### 1. Initial Data Processing
- **`1. initial_data_processing`**: This script handles the initial processing of raw FASTQ files, including merging paired-end reads, dereplicating sequences, and mapping them to the reference genome using Minimap2. It is the first step in the workflow and prepares the data for subsequent analysis.

### 2. Phylogenetic Tree Construction
- **`2. build_tree_NYC.py`**: Builds the phylogenetic tree for the NYC dataset using IQ-TREE.
- **`2. build_tree_ONT.py`**: Constructs the phylogenetic tree for the Ontario dataset using FastTree.

### 3. Extracting Cryptic Branches
- **`3. extract_cryptic_branch_NYC.py`**: Extracts sequences from cryptic branches identified in the NYC phylogenetic tree.
- **`3. extract_cryptic_branch_ONT.py`**: Extracts sequences from cryptic branches identified in the Ontario phylogenetic tree.

### 4. Generating Consensus Sequences
- **`4. generate_consensus_seq_NYC.py`**: Generates a consensus sequence of similar sequences found on the NYC phylogenetic tree.
- **`4. generate_consensus_seq_ONT.py`**: Generates a consensus sequence of similar sequences found on the Ontario phylogenetic tree.

### Common Scripts
- **`derep.py`**: Used in the initial data processing script to dereplicate sequences, ensuring that redundant sequences are removed before further analysis.

## Usage

1. **Start with the initial data processing**: Run `1. initial_data_processing` to process your raw FASTQ files. This script will prepare the data for tree construction by merging, dereplicating, and mapping sequences.

2. **Build the phylogenetic tree**: Depending on your dataset (NYC or Ontario), run either `2. build_tree_NYC.py` or `2. build_tree_ONT.py` to construct the phylogenetic tree.

3. **Extract cryptic branches**: After the tree is built, run the appropriate `3. extract_cryptic_branch_NYC.py` or `3. extract_cryptic_branch_ONT.py` script to extract sequences from cryptic branches.

4. **Generate consensus sequences**: Finally, run `4. generate_consensus_seq_NYC.py` or `4. generate_consensus_seq_ONT.py` to create a consensus sequence from similar sequences identified in the tree.

## Installation and Dependencies

Ensure that you have the following tools and libraries installed:
- Python 3.x
- FastTree (for Ontario)
- IQ-TREE (for NYC)
- MAFFT
- Minimap2
- SAMtools
- Biopython
- CD-HIT
- Additional dependencies listed in `requirements.txt`

## Contributing

Contributions to improve or extend the workflows are welcome! Please submit a pull request or open an issue for any bugs or enhancement requests.

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact

For questions or further information, please contact Veedhi Solanki at veedhisolanki13579@gmail.com.

---

This README now includes a clear note at the beginning that highlights the similarity between the files for both workflows, with the exception of the tree-building step.
