import pickle
import re
import pandas as pd
import os

base = '/home/albertotr/OneDrive/Data/Cambridge_Project/'

# Load the dataframe with the replacements

Lucy_replacements = pd.read_excel(f'{base}/Metadata_genomes.xlsx', sheet_name='Lucy_keys', header=0)

# Debugging the DataFrame
print("DataFrame preview:")
print(Lucy_replacements.head())
print("Columns:", Lucy_replacements.columns)

replacements = {}
replacements_2 = {}
replacements_3 = {}
# Create a dictionary with the replacements
# Iterate through the dataframe and create the replacements dictionary
for index, row in Lucy_replacements.iterrows():
    oldname = str(row['trelabels']).replace('Sample_', '')
    amendment = str(row['Lnames']) + '_' + str(row['year.sampling']).replace(' ', '_')
    camille_name = str(row['final name'])
    replacements.update({oldname: amendment})
    replacements_2.update({amendment: camille_name})
    replacements_3.update({oldname: camille_name})

# Debugging the dictionaries
print("Replacements dictionary:", replacements)
print("Replacements_2 dictionary:", replacements_2)
print("Replacements_3 dictionary:", replacements_3)

sorted_replacements = dict(sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True))
sorted_replacements_2 = dict(sorted(replacements_2.items(), key=lambda x: len(x[0]), reverse=True))
sorted_replacements_3 = dict(sorted(replacements_3.items(), key=lambda x: len(x[0]), reverse=True))

# Debugging the sorted dictionaries
print("Sorted replacements:", sorted_replacements)
print("Sorted replacements_2:", sorted_replacements_2)
print("Sorted replacements_3:", sorted_replacements_3)

# Save the replacements dictionary to a pickle file
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements_2, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save the replacements dictionary to a pickle file
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save the replacements dictionary to a pickle file
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements_foldername.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements_3, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Debugging file paths
print("Base path exists:", os.path.exists(base))