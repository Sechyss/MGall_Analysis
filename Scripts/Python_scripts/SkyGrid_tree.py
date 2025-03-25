#%% Load packages
from ete3 import Tree, TreeStyle, NodeStyle
import matplotlib.pyplot as plt
import pandas as pd
import os

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/')

# Load the tree
tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold_50_highburnin.finaltree.nwk')

metadata = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage_metadata.csv', index_col=0)

# Set the times of the leaves
for leaf in tree_file:
    if leaf.name in metadata.index:
        leaf.add_feature("Year", metadata.loc[leaf.name, "Year"])

# Set the tree style
ts = TreeStyle()
ts.show_leaf_name = True

# Render the tree
tree_file.render("time_phylogenetic_tree.png", tree_style=ts)
# %%
