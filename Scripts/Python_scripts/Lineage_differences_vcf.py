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

clusters_df = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/Visual_60thres_all_VA94/Visual_60thres_all_VA94_microreact_clusters.csv')
lineage1 = list(clusters_df[clusters_df['Rank_6_Lineage_Cluster__autocolour'] == 1]['id'])
lineage2 = list(clusters_df[clusters_df['Rank_6_Lineage_Cluster__autocolour'] == 2]['id'])

lineage_1_snps = non_synonymous[non_synonymous['Sample'].isin(lineage1)]
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

fastafile_creation('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/Different_lineages_VA94_Genes_snps_l1.faa',
genes_affected_l1, '/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCA_000286675.1/genomic.gbff')

#%% COG analysis

#Load the COG dictionary and the hits
cog_dict = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/lineage_differences/'
                         'MM_nma4zs2c.emapper.annotations.xlsx', engine='openpyxl', sheet_name='Sheet2', header=None)
cog_hits = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/lineage_differences/'
                         'MM_nma4zs2c.emapper.annotations.xlsx', engine='openpyxl', sheet_name='Sheet1')

# Create a dictionary from the first and second columns
cog_dictionary = dict(zip(cog_dict[0], cog_dict[1]))

cog_categories = cog_hits['COG_category'].to_list()

data_histogram = {'Function Unknown': 0}
for item in cog_categories:
    if item in cog_dictionary.keys():
        x_values = cog_dictionary[item]
        if x_values not in data_histogram.keys():
            data_histogram[x_values] = 1
        else:
            data_histogram[x_values] += 1

    else:
        data_histogram['Function Unknown'] += 1

# Convert the dictionary to a DataFrame
data_df = pd.DataFrame(list(data_histogram.items()), columns=['COG_function', 'Count'])
# Ensure the COG_function column is of type string
data_df['COG_function'] = data_df['COG_function'].astype(str)

#%% Plotting
# Define a custom color palette
custom_palette = sns.color_palette("husl", len(data_df))
plt.figure(figsize=(10, 6))  # Adjust the figure size

# Create the bar plot with the custom color palette
sns.barplot(x='COG_function', y='Count', hue='COG_function', data=data_df, palette=custom_palette, legend=False)
plt.xticks(fontsize=6, fontweight='bold', rotation='vertical')
plt.yticks(fontweight='bold')
# Rotate the x-axis labels
#plt.xticks(rotation=45)

# Add titles for the x and y axes
plt.xlabel('COG Function')
plt.ylabel('Count')
# Adjust the bottom margin
plt.gcf().subplots_adjust(bottom=0.3)
plt.tight_layout()
plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/COG_function_lineage1.png', dpi=600)
plt.show()