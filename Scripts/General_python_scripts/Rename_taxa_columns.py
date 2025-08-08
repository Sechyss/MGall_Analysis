import os
import pickle
import argparse
import pandas as pd

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

    # Read the Rtab file into a DataFrame
    try:
        df = pd.read_csv(params.file, delimiter='\t')

        # Strip any leading/trailing spaces from the column names
        df.columns = df.columns.str.strip()

        # Replace each key in the dictionary with its corresponding value in the column names
        df.rename(columns=sorted_replacements, inplace=True)

        # Write the modified DataFrame to a new Rtab file
        df.to_csv(new_filename, sep='\t', index=False)
    except FileNotFoundError:
        print(f"Error: The file {params.file} does not exist.")
    except PermissionError:
        print(f"Error: You do not have permission to read/write the file {params.file}.")
