#%% import of packages and function
import pandas as pd


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


#%%

key_data = pd.read_excel('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Metadata_genomes.xlsx',
                         sheet_name='Metadata_Keys')

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

Presence_absence = pd.read_csv('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                               'pangenome_results_filtered/gene_presence_absence.csv')

#%% Filter the relevant sample columns for house finch and poultry

filtered_genes, filtered_genes_2 = filter_presence_absence(Presence_absence, housefinch_all, poultry_all, 0.9, 0.1)

#%% Testing pre and post 2007

filtered_genes_pre2007, filtered_genes_post2007 = filter_presence_absence(Presence_absence, housefinch_pre2007,
                                                                          housefinch_post2007, 0.4, 0.1)
