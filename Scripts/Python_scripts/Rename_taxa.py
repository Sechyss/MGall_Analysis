import os
import pandas as pd

# Define a dictionary with the replacements
dataset_keys = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Metadata_genomes.xlsx', sheet_name='Lucy_keys')

replacements = {}
for index, row in dataset_keys.iterrows():
    oldname = str(row['trelabels']).replace('Sample_', '')
    amendment = str(row['Lnames']) + '_' + str(row['year.sampling'])
    replacements.update({oldname: amendment})

sorted_replacements = dict(sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True))

# Change directory to where the .contree files are located
os.chdir('/home/albertotr/OneDrive/Data/Ipoutcha_motility/BLAST_db/Phylogenetic_trees/')

# Iterate over all files ending with .contree in the directory
for filename in os.listdir():
    if filename.endswith('.contree'):
        # Construct the new filename with 'Edited_' prefix
        new_filename = 'Edited_' + filename

        # Open the current file and create a new file with the edited name
        with open(filename, 'r') as infile, open(new_filename, 'w') as outfile:
            # Iterate through each line in the input file
            for line in infile:
                # Replace each key in the dictionary with its corresponding value
                for old_string, new_string in sorted_replacements.items():
                    line = line.replace(old_string, new_string)
                # Write the modified line to the output file
                outfile.write(line)
