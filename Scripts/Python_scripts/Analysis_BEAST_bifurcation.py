import pandas as pd
import pickle

# Load the CSV file into a DataFrame
Presence_absence = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                               'pangenome_results_filtered/gene_presence_absence.csv')
gene_data = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                               'pangenome_results_filtered/gene_data.csv')

# Load the dictionary of replacement keys from an Excel file
dataset_keys = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Metadata_genomes.xlsx', sheet_name='Lucy_keys')

# Create a dictionary for renaming columns
replacements = {}
for index, row in dataset_keys.iterrows():
    oldname = str(row['trelabels']).replace('Sample_', '')
    amendment = str(row['Lnames']) + '_' + str(row['year.sampling'])
    replacements.update({oldname: amendment})

# Sort replacements by length of key (longest first) to avoid partial replacements
sorted_replacements = dict(sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True))

# Save the sorted replacements dictionary
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/rename_taxa_dict.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Rename the columns in the DataFrame using the replacement dictionary
def replace_column_names(col_name, replacements_dict):
    # Iterate over each key-value pair in the sorted dictionary
    for old, new in replacements_dict.items():
        if old in col_name:
            col_name = col_name.replace(old, new)  # Replace old string with new string
    return col_name

# Apply the renaming function to the column names of the DataFrame
Presence_absence.columns = [replace_column_names(col, sorted_replacements) for col in Presence_absence.columns]

# Check the bifurcation branches
extinct_branch = ['A021_2011', 'A006_2011', 'AMGO_2011', 'A013_2011',
                  'A018_2011', 'A044_2011', 'A004_2011', 'A001_2011', 'A012_2011', 'A3316_2012', 'Blue_2012']
current_branch = ['G_2015', 'L_2015', 'J_2015', 'A020_2011', 'WU47_2013',
                  'OY79_2013', 'A_2015', 'Q_2015', 'Black_2012', 'GW77_2013', 'A1_2014', 'KB64_2013', 'E054_2013',
                  'A046_2011', 'A471_200106_2002244_2012', 'A3278_2012', 'A3240_2012', 'A3225_2012', 'A2486_2012']

genes_old_branch = Presence_absence[extinct_branch].dropna(how='all')
genes_new_branch = Presence_absence[current_branch].dropna(how='all')

oldbranch_unique_genes = set(genes_old_branch.index).difference(set(genes_new_branch.index))
newbranch_unique_genes = set(genes_new_branch.index).difference(set(genes_old_branch.index))

Old_branch_genes = Presence_absence.loc[list(oldbranch_unique_genes)]['Annotation']
New_branch_genes = Presence_absence.loc[list(newbranch_unique_genes)]['Annotation']

