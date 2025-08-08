import os
import pandas as pd
import pickle
import argparse

if '__main__' == __name__:

    parser = argparse.ArgumentParser()

    parser.add_argument('--directory', required=True, type=str)
    params = parser.parse_args()

    # Define a dictionary with the replacements
    dataset_keys = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Metadata_genomes.xlsx', sheet_name='Lucy_keys')

    replacements = {}
    for index, row in dataset_keys.iterrows():
        oldname = str(row['trelabels']).replace('Sample_', '')
        amendment = str(row['Lnames']) + '_' + str(row['year.sampling']).replace(' ', '_')
        replacements.update({oldname: amendment})

    sorted_replacements = dict(sorted(replacements.items(), key=lambda x: len(x[0]), reverse=True))

    with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'wb') as handle:
        pickle.dump(sorted_replacements, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Change directory to where the .contree files are located
    try:
        os.chdir(params.directory)
    except FileNotFoundError:
        print(f"Error: The directory {params.directory} does not exist.")
        exit(1)
    except PermissionError:
        print(f"Error: You do not have permission to access the directory {params.directory}.")
        exit(1)

    for filename in os.listdir():
        if filename.endswith('SRA.fasta'):
            # Construct the new filename with 'Edited_' prefix
            new_filename = 'Edited_' + filename

            # Open the current file and create a new file with the edited name
            try:
                with open(filename, 'r') as infile, open(new_filename, 'w') as outfile:
                    # Iterate through each line in the input file
                    for line in infile:
                        # Replace each key in the dictionary with its corresponding value
                        for old_string, new_string in sorted_replacements.items():
                            if old_string in line:
                                line = line.replace(old_string, new_string)
                                break  # Stop searching for more keys after the first match
                        # Write the modified line to the output file
                        outfile.write(line)
            except FileNotFoundError:
                print(f"Error: The file {filename} does not exist.")
            except PermissionError:
                print(f"Error: You do not have permission to read/write the file {filename}.")
