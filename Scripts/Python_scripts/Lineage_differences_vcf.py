#%%
import os

import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#%%
os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/Snp_differences_lineages/')

snps_lucy = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/All_mutations_matrix.csv', index_col=0)
snps_sra = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/All_mutation_matrix.xlsx', index_col=0, sheet_name='Sheet1', engine='openpyxl')

lucy_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'rb'))
sra_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/SRA_replacements.pickle', 'rb'))

snps_lucy['Sample'] = snps_lucy['Sample'].replace(lucy_replacement)
snps_sra['Sample'] = snps_sra['Sample'].replace(sra_replacement)

all_snps = pd.concat([snps_lucy, snps_sra], axis=0)

non_synonymous = all_snps[all_snps['Effect'].str.contains('STOP') |(all_snps['Effect'].str.contains('NON_SYNONYMOUS')) | (all_snps['Effect'].str.contains('LOST'))]


#%% Lineages

lineage2 = ['S11_1994',
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
'MG30_AL_11_2011'
            ]

lineage_1_snps = non_synonymous[~non_synonymous['Sample'].isin(lineage2)]
lineage_2_snps = non_synonymous[non_synonymous['Sample'].isin(lineage2)]

unique_to_lineage1 = lineage_1_snps[~lineage_1_snps['Position'].isin(lineage_2_snps['Position'])]
unique_to_lineage2 = lineage_2_snps[~lineage_2_snps['Position'].isin(lineage_1_snps['Position'])]

#%% Genes

def fastafile_creation(fasta_file, list_genes, reference_file):
    with open(fasta_file, 'a') as handle:
        for genome in SeqIO.parse(reference_file, 'genbank'):
            for feature in genome.features:
                if feature.type == "CDS":  # Find CDS to collect the information
                    locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '_')
                    if locustag in list_genes:
                        if "translation" in feature.qualifiers.keys():
                            aaseq = Seq(feature.qualifiers["translation"][0])
                            fasta_aa = SeqRecord(aaseq, locustag, description='')
                            SeqIO.write(fasta_aa, handle, 'fasta')


genes_affected_l1 = set(lineage_1_snps['Gene_ID'].to_list())
genes_affected_l2 = set(lineage_2_snps['Gene_ID'].to_list())

fastafile_creation('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/Different_lineages_VA94_Genes_snps_l2.faa',
genes_affected_l2, '/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCA_000286675.1/genomic.gbff')

#%% COG analysis

# Remove all the top part of the excel files
#Load the COG dictionary and the hits
cog_dict = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/lineage_differences/'
                         'COG_dict.xlsx', engine='openpyxl', sheet_name='Sheet1', header=None)
cog_hits_1 = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/lineage_differences/'
                         'Lineage1_EGGNOG-mapped.xlsx', engine='openpyxl', sheet_name='Sheet1')
cog_hits_2 = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/lineage_differences/'
                         'Lineage2_EGGNOG-mapped.xlsx', engine='openpyxl', sheet_name='Sheet1')

# Create a dictionary from the first and second columns
cog_dictionary = dict(zip(cog_dict[0], cog_dict[1]))

cog_categories_1 = cog_hits_1['COG_category'].to_list()

data_histogram = {'Function Unknown': 0}
for item in cog_categories_1:
    if item in cog_dictionary.keys():
        x_values = cog_dictionary[item]
        if x_values not in data_histogram.keys():
            data_histogram[x_values] = 1
        else:
            data_histogram[x_values] += 1

    else:
        data_histogram['Function Unknown'] += 1

cog_categories_2 = cog_hits_2['COG_category'].to_list()

data_histogram_2 = {'Function Unknown': 0}
for item in cog_categories_2:
    if item in cog_dictionary.keys():
        x_values = cog_dictionary[item]
        if x_values not in data_histogram_2.keys():
            data_histogram_2[x_values] = 1
        else:
            data_histogram_2[x_values] += 1
    else:
        data_histogram_2['Function Unknown'] += 1


# Create a set of all unique keys
all_keys = set(data_histogram.keys()).union(set(data_histogram_2.keys()))
# Initialize a new DataFrame with these keys
combined_df = pd.DataFrame(all_keys, columns=['COG_function'])

# Populate the DataFrame with counts from both dictionaries
combined_df['Lineage1'] = combined_df['COG_function'].map(data_histogram).fillna(0).astype(int)
combined_df['Lineage2'] = combined_df['COG_function'].map(data_histogram_2).fillna(0).astype(int)

# Normalize the counts to percentages
total_count_1 = combined_df['Lineage1'].sum()
total_count_2 = combined_df['Lineage2'].sum()

combined_df['Lineage1'] = (combined_df['Lineage1'] / total_count_1) * 100
combined_df['Lineage2'] = (combined_df['Lineage2'] / total_count_2) * 100

# Ensure the COG_function column is of type string
combined_df['COG_function'] = combined_df['COG_function'].astype(str)


#%% Plotting
# Melt the DataFrame to long format for seaborn
melted_df = combined_df.melt(id_vars='COG_function', value_vars=['Lineage1', 'Lineage2'], var_name='Count_Type', value_name='Count')

# Define a custom color palette for each group
palette = {
    'Lineage1': 'skyblue',
    'Lineage2': 'lightgreen'
}

# Create the bar plot with the custom color palette
plt.figure(figsize=(10, 6))  # Adjust the figure size
sns.barplot(x='COG_function', y='Count', hue='Count_Type', data=melted_df, palette=palette)
plt.xticks(fontsize=6, fontweight='bold', rotation='vertical')
plt.yticks(fontweight='bold')

# Add titles for the x and y axes
plt.xlabel('COG Function')
plt.ylabel('Relative Frequency (%)')

# Adjust the bottom margin
plt.gcf().subplots_adjust(bottom=0.3)
plt.tight_layout()
plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/COG_function_combined.png', dpi=600)
plt.show()
# %%
