import pickle
import re
import pandas as pd


base = '/home/albertotr/OneDrive/Data/Cambridge_Project/'

# Load the dataframe with the replacements

Lucy_replacements = pd.read_excel(f'{base}/Metadata_genomes.xlsx', sheet_name='Lucy_keys', header=0)

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

sorted_replacements = dict(sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True))
sorted_replacements_2 = dict(sorted(replacements_2.items(), key=lambda x: len(x[0]), reverse=True))
sorted_replacements_3 = dict(sorted(replacements_3.items(), key=lambda x: len(x[0]), reverse=True))

# Save the replacements dictionary to a pickle file
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements_2, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save the replacements dictionary to a pickle file
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements, handle, protocol=pickle.HIGHEST_PROTOCOL)

# Save the replacements dictionary to a pickle file
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements_foldername.pickle', 'wb') as handle:
    pickle.dump(sorted_replacements_3, handle, protocol=pickle.HIGHEST_PROTOCOL)