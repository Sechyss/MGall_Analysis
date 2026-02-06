import os
import pandas as pd
from ete3 import Tree
os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/Brownian_Alberto/')

def trim_metadata(metadata_file, output_file, rows_to_keep):
    # Read the metadata file
    metadata = pd.read_csv(metadata_file, sep='\t')
    # Keep rows based on values in final.name column
    trimmed_metadata = metadata[metadata['final.name'].isin(rows_to_keep)]
    # Save the trimmed metadata to a new file
    trimmed_metadata.to_csv(output_file, sep='\t', index=False)

def read_tree_and_get_leaves(tree_file):
    tree = Tree(tree_file)
    leaves = tree.get_leaves()
    return [leaf.name for leaf in leaves]

def leafs_to_text_file(leaves, output_file):
    with open(output_file, 'w') as f:
        for leaf in leaves:
            f.write(f"{leaf}\n")

def main():
    metadata_file = 'meta.txt'
    output_file = 'meta_L1.txt'
    tree_file = 'Edited_VA94_consensus_L1_trimmed_60threshold_50_combined.finaltree.txt'
    
    leaves = read_tree_and_get_leaves(tree_file)
    trim_metadata(metadata_file, output_file, leaves)
    leafs_to_text_file(leaves, 'leaves_L1.txt')
if __name__ == "__main__":
    main()