#%% Import libraries and load data
"""
Annotate a phylogenetic tree with SNP-distribution pie charts at each node.

Loads per-sample SNP matrices (Lucy and SRA datasets), maps them onto a
Gubbins-filtered phylogenetic tree, and renders pie-chart glyphs at each
internal node and leaf to summarise the proportion of SNPs per functional
category.  The final annotated tree is exported as a PNG image.

Usage:
    Update the file paths at the top of the script (snps_lucy, snps_sra,
    tree file), then run:
        python Phylogeny_snps.py

Outputs:
    - Annotated tree PNG with per-node SNP category pie charts
"""

from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from collections import Counter
import io
from PIL import Image
import numpy as np

# Load data
snps_lucy = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/All_mutations_matrix.csv', index_col=0)
snps_sra = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/Potential_wrong_run/All_mutation_matrix.xlsx', index_col=0, sheet_name='Sheet1', engine='openpyxl')

lucy_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'rb'))
sra_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/SRA_replacements.pickle', 'rb'))

snps_lucy['Sample'] = snps_lucy['Sample'].replace(lucy_replacement)
snps_sra['Sample'] = snps_sra['Sample'].replace(sra_replacement)

all_snps = pd.concat([snps_lucy, snps_sra], axis=0)

non_synonymous = all_snps[all_snps['Effect'].str.contains('STOP') |(all_snps['Effect'].str.contains('NON_SYNONYMOUS')) | (all_snps['Effect'].str.contains('LOST'))]

mutations = non_synonymous['Effect'].str.split('|').str[3]
non_synonymous['Mutation'] = mutations

#%% Load tree
tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold.finaltree.newick')

# Function to create pie chart
def create_pie_chart(mutations):
    counts = Counter(mutations)
    labels, sizes = zip(*counts.items())
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, autopct='%1.1f%%')
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close(fig)
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, autopct='%1.1f%%')
    buf = io.BytesIO()

# Add pie charts to nodes
for node in tree_file.traverse():
    if not node.is_leaf():
        # Get taxa for this node
        taxa = node.get_leaf_names()
        # Get mutations for these taxa
        node_mutations = non_synonymous[non_synonymous['Sample'].isin(taxa)]['Mutation']
        if not node_mutations.empty:
            pie_chart = create_pie_chart(node_mutations)
            pie_face = faces.ImgFace(pie_chart, width=100, height=100)
            node.add_face(pie_face, column=0, position="branch-right")
        else:
            print(f"No non-synonymous mutations found for taxa: {taxa}")

# Render tree
ts = TreeStyle()
ts.show_leaf_name = True
print(f"No non-synonymous mutations found for taxa: {taxa}")

tree_file.show(tree_style=ts)

# Save the tree to a file
tree_file.render('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_tree_with_pie_charts.png', tree_style=ts)

