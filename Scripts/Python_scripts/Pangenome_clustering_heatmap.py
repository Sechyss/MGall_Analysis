#%% Import libraries and data
"""
Created on Tue Oct 31 14:23:00 2023

This script generates a heatmap of the presence-absence binary matrix of the pangenome and columns are group by hyerarchical clustering.

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

#%% Create the binary matrix
binary_matrix = pd.read_csv(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix.csv', index_col=0)
row_means = binary_matrix.mean(axis=1)
binary_matrix = binary_matrix.loc[~(row_means > 0.95)]

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
# Colour gene clusters for Cas9 protein, Motility genes, and Lipoproteins
cluster_colors = {
    "Cas9 protein": "#d62728",
    "Motility genes": "#2ca02c",
    "Lipoprotein": "#094782",
    "non_cluster": "#7b49b5"
}
row_colors = binary_matrix.index.to_series().map(
    lambda x: (
        cluster_colors["Cas9 protein"] if x in cluster_dict.get('Cas9', []) else
        cluster_colors["Motility genes"] if x in cluster_dict.get('Motility', []) else
        cluster_colors["Lipoprotein"] if x in cluster_dict.get('Lipoprotein', []) else
        cluster_colors["non_cluster"]
    )
)

#%% Create the heatmap and clustering of data
fig, ax = plt.subplots(figsize=(60, 60))
clustering_dend = sns.clustermap(binary_matrix, metric="euclidean", method="ward", 
                                 col_colors=col_colors, row_colors=row_colors,
                                 xticklabels=True, yticklabels=True,
                                 cmap="YlGnBu", 
                                 figsize=(60, 60))
clustering_dend.cax.set_visible(False)
# Reduce the font size of the x-axis labels
clustering_dend.ax_heatmap.set_xticklabels(
    clustering_dend.ax_heatmap.get_xticklabels(), 
    fontsize=10, 
    fontweight='bold'  # Set columns in bold
)
clustering_dend.ax_heatmap.set_yticklabels(
    clustering_dend.ax_heatmap.get_yticklabels(), 
    fontsize=8, 
    fontweight='bold'  # Set rows in bold
)
plt.savefig(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix_clustering.png', dpi=600, bbox_inches='tight')
plt.savefig(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix_clustering.pdf', dpi=600, bbox_inches='tight')
plt.show()
# %%
