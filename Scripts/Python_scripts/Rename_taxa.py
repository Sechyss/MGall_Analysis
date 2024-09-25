import os
import pandas as pd

# Define a dictionary with the replacements
dataset_keys = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Metadata_genomes.xlsx', sheet_name='Lucy_keys')

replacements ={}
for index, row in dataset_keys.iterrows():
    oldname = str(row['trelabels']).replace('Sample_', '')
    amendment = str(row['Lnames'])+'_'+str(row['year.sampling'])
    replacements.update({oldname:amendment})

sorted_replacements = dict(sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True))
os.chdir('/home/albertotr/downloads/')
# Open the input text file in read mode and a new output file in write mode
with open('Rlow_gubbins.masked_Combined.txt', 'r') as infile, open('Rlow_lucy_masked_gubbins_edited.txt', 'w') as outfile:
    # Iterate through each line in the input file
    for line in infile:
        # Replace each key in the dictionary with its corresponding value
        for old_string, new_string in sorted_replacements.items():
            line = line.replace(old_string, new_string)
        # Write the modified line to the output file
        outfile.write(line)
