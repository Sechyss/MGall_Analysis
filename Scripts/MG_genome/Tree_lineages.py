#%% Load phylo_env to run this code, it has ete3 installed
"""
Cluster tree leaves into lineages using patristic distances.

Loads an IQ-TREE consensus tree, computes the full pairwise patristic
distance matrix, applies hierarchical clustering (Ward linkage), and
assigns each leaf to a lineage cluster.  The resulting assignments are
visualised and can be exported for use in downstream analyses.

Requirements:
    Run in an environment that includes ete3 (e.g. activate phylo_env).

Usage:
    Update the tree file path, then run:
        python Tree_lineages.py

Outputs:
    - Dendrogram and cluster assignment plots (displayed or saved as PNG)
    - Lineage-membership dictionary (in-memory or serialised to pickle)
"""

import random

import matplotlib.pyplot as plt
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram

# Load the tree
tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/IQtree_informativesites/Rlowsnps_stripped_informativesites.contree')
tree_file.unroot()
# Get all leaves in the tree
leaves = tree_file.get_leaves()
n = len(leaves)

# Initialize an empty distance matrix
distance_matrix = np.zeros((n, n))

# Calculate pairwise distances between all leaves
for i, leaf1 in enumerate(leaves):
    for j, leaf2 in enumerate(leaves):
        if i < j:  # Only compute upper triangle (distance matrix is symmetric)
            dist = tree_file.get_distance(leaf1, leaf2)
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist  # Mirror the value

# Convert the distance matrix to a condensed form for clustering
condensed_matrix = distance_matrix[np.triu_indices(n, k=1)]

# Perform hierarchical clustering
Z = linkage(condensed_matrix, 'ward')  # You can also use other linkage methods
threshold = 0.6
clusters = fcluster(Z, threshold, criterion="distance")

# Plot the dendrogram with the threshold line
plt.figure(figsize=(10, 8))
dendrogram(Z, labels=[leaf.name for leaf in leaves], leaf_rotation=90)
plt.axhline(y=threshold, color='r', linestyle='--', label=f'Threshold = {threshold}')
plt.legend()
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("Leaf Index or Clustered Sample")
plt.ylabel("Genetic Distance")
plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Hierarchical_lineage/IQTree_snps_Lucy_infosites_Clustering_dendrogram.pdf')

# Display lineage clusters
print(f"Estimated lineages (clusters): {np.unique(clusters)}")

#%% Assign each leaf and branch a cluster color
cluster_colors = {}
for cluster_id in np.unique(clusters):
    # Assign a random color for each cluster
    cluster_colors[cluster_id] = "#" + ''.join([random.choice('0123456789ABCDEF') for _ in range(6)])

# Function to color branches from each leaf up to the common ancestor
def color_branch_upwards(node, color):
    # Apply color to the node's branch
    style = NodeStyle()
    style["fgcolor"] = color          # Set the font color for the leaf name
    style["vt_line_color"] = color    # Set the color for the vertical branch
    style["hz_line_color"] = color    # Set the color for the horizontal branch
    style["size"] = 0                 # Hide the node circle
    node.set_style(style)

    # Continue coloring up to the root
    while node.up:
        node = node.up
        node.set_style(style)

# Apply colors to leaves and their branches based on clusters
for i, leaf in enumerate(leaves):
    cluster_id = clusters[i]
    color = cluster_colors[cluster_id]
    color_branch_upwards(leaf, color)  # Color the leaf and its ancestors
    print(f"Leaf {leaf.name} is in lineage cluster {cluster_id} with color {color}")

# Define a TreeStyle for display
tree_style = TreeStyle()
tree_style.show_leaf_name = True
tree_style.scale = 20  # Adjust scale for visibility
tree_style.branch_vertical_margin = 10  # Increase spacing between branches

# Render and save the tree to a PNG file
tree_file.show(tree_style=tree_style)  # Increased width and height for better resolution
