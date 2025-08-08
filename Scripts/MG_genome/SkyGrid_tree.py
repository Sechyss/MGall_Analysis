#%% Load packages
from ete3 import Tree, TreeStyle, NodeStyle, Nexml
import pandas as pd
import os

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/')

# Load the tree
tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold_50_highburnin.finaltree.nwk')
nexml_project = Nexml()
#Load nexml file
nexml_project.build_from_file('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold_50_highburnin.finaltree.nexus')
t = nexml_project.trees()[0]
metadata = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage_metadata.csv', index_col=0)

# Set the times of the leaves
for leaf in tree_file:
    if leaf.name in metadata.index:
        leaf.add_feature("Year", metadata.loc[leaf.name, "Year"])

# Set the tree style
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = True
ts.show_branch_support = True
ts.show_scale = True

# Set the node style for group 1
nstyle = NodeStyle()
nstyle['hz_line_color'] = '#56adf5'
nstyle['vt_line_color'] = '#56adf5'
nstyle['hz_line_width'] = 2
nstyle['vt_line_width'] = 2
for node in tree_file.traverse():
    node.set_style(nstyle)
#nstyle['hz_line_color'] = '#f5ad56'

tree_file.show(tree_style=ts)

# Render the tree
#tree_file.render("time_phylogenetic_tree.png", tree_style=ts, units="px", w=1000, dpi=600, h=1000)
# %%
