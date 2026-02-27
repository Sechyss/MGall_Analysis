#!/usr/bin/env python3
"""
Replace taxon names in any text-based file using a serialised dictionary.

Loads a pickle dictionary that maps old sample/taxon identifiers to new ones
and performs whole-word substitutions throughout a text file (e.g. Newick
trees, alignment headers, TSV tables).  The edited content is written to a
new file prefixed with 'Edited_'.

Usage:
    python Rename_taxa_file.py --file <input_file> --dictionary <replacements.pickle>

Arguments:
    --file        Path to the text file whose taxon names should be replaced.
    --dictionary  Path to a pickle file containing a dict {old_name: new_name}.

Outputs:
    - Edited_<input_filename>  Modified file with updated taxon names
"""
import os
import pickle
import argparse
import re

if '__main__' == __name__:

    parser = argparse.ArgumentParser()

    parser.add_argument('--file', required=True, type=str)
    parser.add_argument('--dictionary', required=True, type=str)
    params = parser.parse_args()

    # Load the dictionary with the replacements
    with open(params.dictionary, 'rb') as handle:
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
                    # Use regular expression to match whole words
                    line = re.sub(rf'\b{re.escape(old_string)}\b', new_string, line)
                # Write the modified line to the output file
                outfile.write(line)
    except FileNotFoundError:
        print(f"Error: The file {params.file} does not exist.")
    except PermissionError:
        print(f"Error: You do not have permission to read/write the file {params.file}.")
