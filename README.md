# Cryptic-Lineage-Detection-in-Wastewater-Surveillance-for-COVID-19
This repository provides workflows and scripts designed for the detection of cryptic SARS-CoV-2 lineages in wastewater samples. The repository includes two main workflows tailored for data from Ontario and NYC, with specialized scripts for each region while maintaining a common underlying structure.

Project Structure
Ontario Workflow
extract_seq.py: Extracts sequences from cryptic branches identified in the phylogenetic tree.
fast2.py: Generates the phylogenetic tree using FastTree, optimized for the Ontario dataset.
consensus.py: Generates a consensus sequence from all similar sequences on the phylogenetic tree.
NYC Workflow
build_tree.py: Builds the phylogenetic tree using IQ-TREE, tailored for the NYC dataset.
extract_long_branch.py: Extracts sequences from long branches, identifying potential cryptic lineages.
generate_consensus_seq.py: Generates a consensus sequence from similar sequences found on the tree.
Common Scripts
derep.py: Used in both workflows to dereplicate sequences, ensuring that redundant sequences are removed before further analysis.
Usage
Data Preprocessing: Start with derep.py to dereplicate your sequence data. This step is common for both Ontario and NYC workflows.

Phylogenetic Tree Generation:

For Ontario data, use fast2.py to generate the tree using FastTree.
For NYC data, use build_tree.py to generate the tree using IQ-TREE.
Extracting Cryptic Lineages:

Use extract_seq.py (Ontario) or extract_long_branch.py (NYC) to extract sequences from branches suspected to contain cryptic lineages.
Consensus Sequence Generation:

Run consensus.py (Ontario) or generate_consensus_seq.py (NYC) to produce a consensus sequence of similar sequences on the tree.
Installation and Dependencies
Ensure that you have the following tools and libraries installed:

Python 3.x
FastTree (for Ontario)
IQ-TREE (for NYC)
Biopython
Additional dependencies listed in requirements.txt
Contributing
Contributions to improve or extend the workflows are welcome! Please submit a pull request or open an issue for any bugs or enhancement requests.

License
This project is licensed under the MIT License. See the LICENSE file for more details.

Contact
For questions or further information, please contact Veedhi Solanki at veedhisolanki13579@gmail.com

