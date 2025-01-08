#%% import of packages and function
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pickle


def filter_presence_absence(dataframe, group1, group2, filter1, filter2):
    dataframe_group1 = [col for col in dataframe.columns if col in group1]
    dataframe_group2 = [col for col in dataframe.columns if col in group2]

    # Calculate the presence ratio for house finch and poultry
    presence_group1 = dataframe[dataframe_group1].notna().mean(axis=1)
    presence_group2 = dataframe[dataframe_group2].notna().mean(axis=1)

    # Filter genes present in at least 90% of house finch samples and less than 10% of poultry samples
    filtered_genes_group1 = dataframe[
        (presence_group1 >= filter1) & (presence_group2 <= filter2)]
    filtered_genes_group2 = dataframe[
        (presence_group1 <= filter2) & (presence_group2 >= filter1)]
    return filtered_genes_group1, filtered_genes_group2

def fastafile_creation(fasta_file, list_genes, reference_file):
    pangenome_reference = SeqIO.parse(reference_file, 'fasta')
    with open(fasta_file+'.fasta', 'a') as f2:
        with open(fasta_file+'.faa', 'a') as f1:
            for gene in pangenome_reference:
                if gene.id in list_genes:
                    sequence = gene.seq
                    prot_sequence = sequence.translate(table=4, stop_symbol='')
                    seq_record = SeqRecord(sequence, id=gene.id, description='')
                    seq_record_2 = SeqRecord(prot_sequence, id=gene.id, description='')
                    SeqIO.write(seq_record, f2, 'fasta')
                    SeqIO.write(seq_record_2, f1, 'fasta')

Presence_absence = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                               'pangenome_results_filtered/gene_presence_absence.csv')

#%% Set up the data

key_data = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Metadata_genomes.xlsx',
                         sheet_name='Metadata_Keys')
key_data['Sample Name'] = key_data['Sample Name'].replace(
    {'/': '_', ' ': '_', '\\.': '_', '-': '_', '\(': '_', '\)': ''}
    , regex=True)

housefinch_all = list(key_data[key_data['Host'].isin(['Haemorhous mexicanus', 'House finch'])]['Sample Name'])
poultry_all = list(key_data[~key_data['Sample Name'].isin(housefinch_all)]['Sample Name'])

housefinch_pre2007 = list(
    key_data[(key_data['Host'].isin(['Haemorhous mexicanus', 'House finch'])) & (key_data['Date'] < 2007)][
        'Sample Name'])
housefinch_post2007 = list(
    key_data[(key_data['Host'].isin(['Haemorhous mexicanus', 'House finch'])) & (key_data['Date'] >= 2007)][
        'Sample Name'])

poultry_pre2007 = list(
    key_data[(~key_data['Host'].isin(['Haemorhous mexicanus', 'House finch'])) & (key_data['Date'] < 2007)][
        'Sample Name'])
poultry_post2007 = list(
    key_data[(~key_data['Host'].isin(['Haemorhous mexicanus', 'House finch'])) & (key_data['Date'] >= 2007)][
        'Sample Name'])

filtering_dict = {'Poultry_all': poultry_all, 'Poultry_pre2007': poultry_pre2007, 'Poultry_post2007': poultry_post2007,
                  'Housefinch_all': housefinch_all, 'Housefinch_post2007': housefinch_post2007,
                  'Housefinch_pre2007': housefinch_pre2007}
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/filtering_dict.pickle', 'wb') as handle:
    pickle.dump(filtering_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

#%% Filter the relevant sample columns for house finch and poultry

filtered_genes, filtered_genes_2 = filter_presence_absence(Presence_absence, housefinch_all, poultry_all, 0.9, 0.1)

#%% Testing pre and post 2007

filtered_genes_pre2007, filtered_genes_post2007 = filter_presence_absence(Presence_absence, housefinch_pre2007,
                                                                          housefinch_post2007, 0.15, 0)

#%% Differences between poultry and early house finch and poultry

filtered_genes_early_HF, filtered_genes_poultry = filter_presence_absence(Presence_absence, housefinch_pre2007,
                                                                          poultry_pre2007, 0.15, 0)

#%% Analysis presence/absence

lineages_table = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/mgall_60threshold/'
                             'mgall_60threshold_lineages.csv')
lineages= dict(zip(lineages_table['id'], lineages_table['Rank_2_Lineage']))
lineages_dict = {key.split('_')[0]: value for key, value in lineages.items()}

key_data = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Metadata_genomes.xlsx',
                         sheet_name='Lucy_keys')
translation_dict = dict(zip(key_data['Lnames'], key_data['trelabels']))
translated_dict = {}
for keys in translation_dict.keys():
    newkey = lineages_dict[keys]
    if newkey not in translated_dict.keys():
        translated_dict[newkey] = [translation_dict[keys].replace('Sample_', '')]
    else:
        translated_dict[newkey].append(translation_dict[keys].replace('Sample_', ''))


with open('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/mgall_60threshold/lineage_dict.pickle', 'wb') as handle:
    pickle.dump(translated_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

lineage1 = set(translated_dict[1])
lineage2 = set(translated_dict[2])



#%% Lineages Rlow

lineages_table_Rlow = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/mgall_60threshold/'
                                  'mgall_60threshold_lineages.csv')
lineages_Rlow= dict(zip(lineages_table_Rlow['id'], lineages_table_Rlow['Rank_2_Lineage']))
lineages_dict_Rlow = {key.split('_')[0]: value for key, value in lineages_Rlow.items()}

translated_dict_Rlow = {}
for keys in translation_dict.keys():
    newkey = lineages_dict_Rlow[keys]
    if newkey not in translated_dict_Rlow.keys():
        translated_dict_Rlow[newkey] = [translation_dict[keys].replace('Sample_', '')]
    else:
        translated_dict_Rlow[newkey].append(translation_dict[keys].replace('Sample_', ''))

lineage1_Rlow = set(translated_dict_Rlow[1])
lineage2_Rlow = set(translated_dict_Rlow[2])
lineage3_Rlow = set(translated_dict_Rlow[3])

lineage2_Rlow_relative = lineage2_Rlow.union(lineage3_Rlow)

#%% Venn diagram of lineages
from venny4py.venny4py import *

sets = {
    'Lineage1_Rlow': lineage1_Rlow,
    'Lineage2_Rlow': lineage2_Rlow_relative,
    'Lineage1_VA94': lineage1,
    'Lineage2_VA94': lineage2
}

venny4py(sets=sets)

plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/Lineages_60threshold_venn.png', dpi=600)

#%% filtering by lineage

lineage1_genes, lineage2_genes = filter_presence_absence(Presence_absence, lineage1, lineage2, 0.15, 0.0)

list_of_genes = lineage1_genes['Gene'].to_list()

pangenome_file = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/pan_genome_reference.fa'
output_fasta = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/Lineage_1_shell_differences'
fastafile_creation(output_fasta, list_of_genes, pangenome_file)