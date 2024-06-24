import pandas as pd

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
                               'pangenome_results_genomes/gene_presence_absence.csv')

# Filter the relevant sample columns for house finch and poultry
housefinch_samples = [col for col in Presence_absence.columns if col in housefinch_all]
poultry_samples = [col for col in Presence_absence.columns if col in poultry_all]

# Calculate the presence ratio for house finch and poultry
housefinch_presence = Presence_absence[housefinch_samples].notna().mean(axis=1)
poultry_presence = Presence_absence[poultry_samples].notna().mean(axis=1)

# Filter genes present in at least 90% of house finch samples and less than 10% of poultry samples
filtered_genes = Presence_absence[(housefinch_presence >= 0.9) & (poultry_presence < 0.1)]
filtered_genes_2 = Presence_absence[(housefinch_presence < 0.1) & (poultry_presence >= 0.9)]
