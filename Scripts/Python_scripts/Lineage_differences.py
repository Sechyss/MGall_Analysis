#%% Importing the libraries
import pickle
import pandas as pd
import math

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from venny4py.venny4py import *
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from ete3 import Tree
import logging

logging.basicConfig(level=logging.INFO)

def filter_presence_absence(dataframe, group1, group2, filter1, filter2):
    dataframe_group1 = [col for col in dataframe.columns if col in group1]
    dataframe_group2 = [col for col in dataframe.columns if col in group2]

    # Calculate the presence ratio for group1 and group2
    presence_group1 = dataframe[dataframe_group1].notna().mean(axis=1)
    presence_group2 = dataframe[dataframe_group2].notna().mean(axis=1)

    # Filter genes present in at least filter1 threshold of group1 samples and less than filter2 threshold of group2 samples
    filtered_genes_group1 = dataframe[
        (presence_group1 >= filter1) & (presence_group2 < filter2)]
    filtered_genes_group2 = dataframe[
        (presence_group1 < filter2) & (presence_group2 >= filter1)]
    return filtered_genes_group1, filtered_genes_group2


base_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_HF'
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

lucy_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'rb'))
sra_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/SRA_replacements.pickle', 'rb'))
#%% Load the data and the lineages
snps_lucy = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/All_mutations_matrix.csv', index_col=0)
snps_sra = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/All_mutation_matrix.xlsx', index_col=0, sheet_name='Sheet1', engine='openpyxl')

snps_lucy['Sample'] = snps_lucy['Sample'].replace(lucy_replacement)
snps_sra['Sample'] = snps_sra['Sample'].replace(sra_replacement)

all_snps = pd.concat([snps_lucy, snps_sra], axis=0)

non_synonymous = all_snps[all_snps['Effect'].str.contains('STOP') |(all_snps['Effect'].str.contains('NON_SYNONYMOUS')) | (all_snps['Effect'].str.contains('LOST'))]

lineage_1_snps = all_snps[all_snps['Sample'].isin(lineage1)]
lineage_2_snps = all_snps[all_snps['Sample'].isin(lineage2)]

unique_to_lineage1 = lineage_1_snps[~lineage_1_snps['Position'].isin(lineage_2_snps['Position'])]
unique_to_lineage2 = lineage_2_snps[~lineage_2_snps['Position'].isin(lineage_1_snps['Position'])]

genes_affected_l1 = set(lineage_1_snps['Gene_ID'].to_list())
genes_affected_l2 = set(lineage_2_snps['Gene_ID'].to_list())

unique_genes_lineage2 = genes_affected_l2.difference(genes_affected_l1)
filtered_df = lineage_2_snps[lineage_2_snps['Gene_ID'].isin(unique_genes_lineage2)]

#%% Venn diagram to show overlapping snps between lineages

sets = {
    'Lineage1_genes': {x for x in genes_affected_l1 if not (isinstance(x, float) and math.isnan(x))},
    'Lineage2_genes': {x for x in genes_affected_l2 if not (isinstance(x, float) and math.isnan(x))},
}

venny4py(sets=sets, out='/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_HF/Lineage_differences/Different_lineages_VA94_Genes_snps_Venn.png',dpi=600)

#%% Compare presence-absence of genes from the pangenome

presence_absence_file = f'{base_path}/gene_presence_absence_filt_pseudo_length_frag.csv'
venn_output_path = f'{base_path}/Lineage_differences/Different_lineages_VA94_Genes_snps_Venn.png'
heatmap_output_path = f'{base_path}/Lineage_differences/gene_presence_absence_heatmap_with_gap.png'

try:
    Presence_absence = pd.read_csv(presence_absence_file)
except FileNotFoundError:
    print(f"File not found: {presence_absence_file}")
    raise

replacements = {**lucy_replacement, **sra_replacement}
Presence_absence.columns = Presence_absence.columns.to_series().replace(replacements)

lineage1_df, lineage2_df = filter_presence_absence(Presence_absence, lineage1, lineage2, 0.60, 0.30)
lineage1_genes = lineage1_df['Gene'].to_list()
lineage2_genes = lineage2_df['Gene'].to_list()
#%% Fasta file creation of the genes in the different lineages

def fastafile_creation(fasta_file, list_genes, reference_file):
    with open(fasta_file, 'a') as handle:
        fasta_db = SeqIO.parse(reference_file, 'fasta')
        for record in fasta_db:
            if record.id in list_genes:
                fasta_aa = SeqRecord(record.seq, record.id, description='')
                SeqIO.write(fasta_aa, handle, 'fasta')
    handle.close()

fastafile_creation(f'{base_path}/Lineage_differences/Lineage1_genes.fna', lineage1_genes, f'{base_path}/pan_genome_reference.fa')
fastafile_creation(f'{base_path}/Lineage_differences/Lineage2_genes.fna', lineage2_genes, f'{base_path}/pan_genome_reference.fa')

#%% Heatmap of the presence-absence of genes in the two lineages
FIG_SIZE = (30, 20)
FONT_SIZE = 20
COLORS = ['#1f77b4', '#ff7f0e']  # Blue and Orange

restructured_df = Presence_absence.set_index('Gene')
# Reorder the columns to have lineage 1 first and then lineage 2
ordered_columns = lineage1 + lineage2
restructured_df = restructured_df[ordered_columns]
# Create a binary matrix (1 for presence, 0 for absence)
binary_matrix = restructured_df.notna().astype(int)
binary_matrix.to_csv(f'{base_path}/Lineage_differences/Presence_absence_binary_matrix.csv')
row_means = binary_matrix.mean(axis=1)
binary_matrix = binary_matrix.loc[~(row_means > 0.95)]

# Create the heatmap with separation between lineages

# Split the binary matrix into two parts for lineage 1 and lineage 2
binary_matrix_lineage1 = binary_matrix[lineage1]
binary_matrix_lineage2 = binary_matrix[lineage2]

# Add a column of NaN values to create a gap between the two lineages
gap = pd.DataFrame(np.nan, index=binary_matrix.index, columns=['Gap'])
binary_matrix_with_gap = pd.concat([binary_matrix_lineage1, gap, binary_matrix_lineage2], axis=1)

col_colors = sns.color_palette(COLORS)  # Blue for lineage 1, Orange for lineage 2

# Plot the heatmap
plt.figure(figsize=FIG_SIZE)  # Adjust the figure size as needed
sns.heatmap(binary_matrix_with_gap, cmap=col_colors, cbar=False, linewidths=0.5, yticklabels=False)

# Update x-axis labels to exclude "Gap"
xtick_labels = lineage1 + [''] + lineage2  # Replace "Gap" with an empty string
plt.xticks(ticks=np.arange(len(xtick_labels)) + 0.5, labels=xtick_labels, rotation=90)

plt.xlabel('Samples', fontsize=FONT_SIZE)
plt.ylabel('Genes', fontsize=FONT_SIZE)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(heatmap_output_path, dpi=600)
plt.show()

logging.info("Heatmap created successfully.")

# %%
