#%% Import libraries and data
"""
Generate a hierarchical clustering heatmap of the pangenome presence/absence matrix.

Reads the Panaroo gene presence/absence matrix, orders samples by lineage
membership, performs hierarchical clustering on both genes and samples, and
produces an annotated heatmap.  A cluster-category dictionary (from
Dictionary_clusters_heatmap.py) is used to colour-annotate gene rows by
functional category (lipoproteins, virulence, motility, etc.).

Usage:
    Update 'base_path' and associated paths, then run:
        python Pangenome_clustering_heatmap.py

Outputs:
    - Heatmap PNG of gene presence/absence, clustered by gene and sample
    - Sorted gene list CSV for downstream filtering
"""
import pandas as pd
import seaborn as sns
from ete3 import Tree
import pickle
from matplotlib import pyplot as plt

# Create the lineages
lineage2 = [
'A090809_2009',
'G_2015',
'L_2015',
'J_2015',
'A01546_2007',
'A01505_2007',
'A01510_2007',
'A565_2002',
'A2261_2000',
'A506_2002',
'A454_2001',
'A442_2001',
'A471_2001',
'A267_2000',
'A3611_2001',
'A509_2002',
'A277_2001',
'BP421_2002',
'A933_2001',
'A297_2001',
'A371_2001',
'A877_2003',
'MG31_AL_11_2011',
'A020_2011',
'OY79_2013',
'WU47_2013',
'A_2015',
'Q_2015',
'Black_2012',
'GW77_2013',
'E054_2013',
'KB64_2013',
'A1_2014',
'MG27_AL_11_2011',
'A046_2011',
'A3278_2012',
'A3240_2012',
'A3244_2012',
'A3225_2012',
'A2486_2012',
'A304_2003',
'A195_2003',
'A012_2011',
'MG28_AL_11_2011',
'A044_2011',
'A001_2011',
'MG22_AL_11_2011',
'A004_2011',
'A3316_2012',
'Blue_2012',
'A021_2011',
'MG23_AL_11_2011',
'A006_2011',
'MG25_AL_11_2011',
'MG32_AL_11_2011',
'AMGO_2011',
'MG29_AL_11_2011',
'A013_2011',
'MG26_AL_11_2011',
'A018_2011',
'MG30_AL_11_2011',
'NC08_2008.031-4-3P',
'NC06_2006.080-5-2P'
            ]

tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold_50_combined.finaltree.newick')
leaves = tree_file.get_leaves()

lineage1 = [leaf.name.replace("'", "") for leaf in leaves if leaf.name.replace("'", "") not in lineage2]
lineage1 = [elem.replace('_2011', '') if 'MG' in elem and '_AL_11_2011' in elem else elem for elem in lineage1]
lineage2 = [elem.replace('_2011', '') if 'MG' in elem and '_AL_11_2011' in elem else elem for elem in lineage2]

base_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_HF'


# Import the dictionary with the lipoproteins
with open(f'{base_path}/candidates_clusters_pangenome.pickle', 'rb') as f:
    cluster_dict = pickle.load(f)
# Import the dictionary with the replacements
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements.pickle', 'rb') as f:
    replacements = pickle.load(f)

# Replace lineage1 and lineage2 names based on the replacements dictionary
lineage1 = [replacements.get(name, name) for name in lineage1]
lineage2 = [replacements.get(name, name) for name in lineage2]

#%% Create the binary matrix
binary_matrix = pd.read_csv(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix.csv', index_col=0)
# Replace column names based on the replacements dictionary
binary_matrix.columns = binary_matrix.columns.map(lambda x: replacements.get(x, x))
row_means = binary_matrix.mean(axis=1)
binary_matrix = binary_matrix.loc[~(row_means > 0.95)]
# Reorganize columns to put lineage1 first and then lineage2
lineage1_columns = [col for col in binary_matrix.columns if col in lineage1]
lineage2_columns = [col for col in binary_matrix.columns if col in lineage2]
binary_matrix = binary_matrix[lineage1_columns + lineage2_columns]

# Assign colors to each lineage
lineage_colors = {
    "lineage1": "#1f77b4",
    "lineage2": "#ff7f0e"
}

# Create a color palette for the lineages

# Map colors to columns (samples) based on their lineage
col_colors = binary_matrix.columns.to_series().map(
    lambda x: lineage_colors["lineage1"] if x in lineage1 else lineage_colors["lineage2"]
)
#%% Alternative to colour all as one single group - Virulence genes
cluster_colors = {
    "Virulence_genes": "#3411bf",
    "non_cluster": "#ff5733"
}

row_colors = binary_matrix.index.to_series().map(
    lambda x: (
        cluster_colors["Virulence_genes"] if x in cluster_dict.get('Virulence', []) else
        cluster_colors["non_cluster"]
    )
)

#%% Alternative to colour different clusters
# Create a color palette for the clusters
# Colour gene clusters for Cas9 protein, Virulence related (Motility + Lipoproteins), and others
cluster_colors = {
    "Cas9 protein": "#d62728",
    "Virulence related": "#2ca02c",  # Combined Motility + Lipoprotein
    "non_cluster": "#7b49b5"
}

row_colors = binary_matrix.index.to_series().map(
    lambda x: (
        cluster_colors["Cas9 protein"] if x in cluster_dict.get('Cas9', []) else
        cluster_colors["Virulence related"] if x in cluster_dict.get('Motility', []) or x in cluster_dict.get('Lipoprotein', []) else
        cluster_colors["non_cluster"]
    )
)

#%% Create the heatmap and clustering of data
# Remove the standalone figure creation - clustermap creates its own
clustering_dend = sns.clustermap(binary_matrix, metric="euclidean", method="ward", 
                                 col_colors=col_colors, row_colors=row_colors,
                                 col_cluster=False, row_cluster=True,
                                 xticklabels=True, yticklabels=True,
                                 cmap="YlGnBu", 
                                 figsize=(60, 60),
                                 dendrogram_ratio=(0.1, 0.05),  # Reduce dendrogram size (row, col)
                                 cbar_pos=None)  # Remove colorbar completely

# Adjust the layout to reduce top space
clustering_dend.ax_row_dendrogram.set_visible(True)
clustering_dend.ax_col_dendrogram.set_visible(False)  # Hide column dendrogram since col_cluster=False

# Reduce the font size of the x-axis labels
clustering_dend.ax_heatmap.set_xticklabels(
    clustering_dend.ax_heatmap.get_xticklabels(), 
    fontsize=10, 
    fontweight='bold'
)
clustering_dend.ax_heatmap.set_yticklabels(
    clustering_dend.ax_heatmap.get_yticklabels(), 
    fontsize=8, 
    fontweight='bold'
)

# Add legend for row colors (Gene clusters) - BIGGER FONTS
from matplotlib.patches import Patch
row_legend_handles = [
    Patch(facecolor=cluster_colors["Cas9 protein"], label="Cas9 protein"),
    Patch(facecolor=cluster_colors["Virulence related"], label="Virulence related"),
    Patch(facecolor=cluster_colors["non_cluster"], label="Other genes")
]

# Add legend for column colors (Groups) - BIGGER FONTS
col_legend_handles = [
    Patch(facecolor=lineage_colors["lineage1"], label="Group 1"),
    Patch(facecolor=lineage_colors["lineage2"], label="Group 2")
]

# Combine both legends and place at top
all_handles = row_legend_handles + col_legend_handles
all_labels = ["Cas9 protein", "Virulence related", "Other genes", "Group 1", "Group 2"]

# Create a single legend at the top of the figure
fig = clustering_dend.fig
fig.legend(
    handles=all_handles,
    labels=all_labels,
    loc='upper center',
    bbox_to_anchor=(0.5, 0.98),
    ncol=5,  # Updated: 5 items instead of 6
    frameon=True,
    fontsize=42,
    title="Gene Categories & Sample Groups",
    title_fontsize=50,
    edgecolor='black',
    fancybox=True,
    shadow=True
)

# Adjust subplot positions to reduce top whitespace and make room for legend
clustering_dend.gs.update(top=0.92, bottom=0.05, left=0.05, right=0.95)

plt.savefig(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix_clustering.png', dpi=600, bbox_inches='tight')
plt.savefig(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix_clustering.pdf', dpi=600, bbox_inches='tight')
plt.savefig(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix_clustering.svg', dpi=600, bbox_inches='tight')
