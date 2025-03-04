import os
import pickle
import argparse

if '__main__' == __name__:

    parser = argparse.ArgumentParser()

    parser.add_argument('--file', required=True, type=str)
    params = parser.parse_args()

    # Load the dictionary with the replacements
    with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'rb') as handle:
        sorted_replacements = pickle.load(handle)

    # Check if the file exists
    if not os.path.isfile(params.file):
        print(f"Error: The file {params.file} does not exist.")
        exit(1)

    # Construct the new filename with 'Edited_' prefix
    new_filename = 'Edited_' + os.path.basename(params.file)

    # Open the current file and create a new file with the edited name
    try:
        with open(params.file, 'r') as infile, open(new_filename, 'w') as outfile:
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
        print(f"Error: The file {params.file} does not exist.")
    except PermissionError:
        print(f"Error: You do not have permission to read/write the file {params.file}.")
