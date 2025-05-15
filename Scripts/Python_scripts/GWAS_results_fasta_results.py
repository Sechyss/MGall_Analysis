# usr/bin/env python
# -*- coding: utf-8 -*-
# This script reads GWAS results from a CSV or TSV file, filters the results based on a p-value threshold,
# and extracts the corresponding sequences from a FASTA file. The filtered sequences are then written to a new FASTA file.
# Ensure you have the required libraries installed:
# pip install pandas biopython
# Make sure to replace the file paths with your actual file paths before running the script.
# This script is designed to be run in a Python environment with the necessary libraries installed.
# The script assumes that the GWAS results file contains a column named 'filter-pvalue' for filtering.
# The script also assumes that the FASTA file contains sequence records with IDs matching the variant IDs in the GWAS results.
from matplotlib import table
import pandas as pd
import os
from Bio import SeqIO

def read_table(file_path):
    """
    Read a table from a file and return a DataFrame and extract the cluster information.
    """
    if file_path.endswith('.csv'):
        table = pd.read_csv(file_path, header=0)
        filtered_table = table[table['filter-pvalue'] < 0.05]
        return filtered_table['variant'].tolist()
    elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
        # Assuming tab-separated values for TSV or TXT files
        table = pd.read_csv(file_path, sep='\t', header=0)
        filtered_table = table[table['filter-pvalue'] < 0.05]
        return filtered_table['variant'].tolist()
    else:
        raise ValueError("Unsupported file format. Please provide a .csv nor .tsv.")

def extract_sequences(fasta_file, variant_ids):
    """
    Extract sequences from a FASTA file based on the provided variant IDs.
    """
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in variant_ids:
            sequences.append(record)
    return sequences

def write_fasta(sequences, output_file):
    """
    Write the extracted sequences to a new FASTA file.
    """
    with open(output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
def main():
    # Define the input files
    os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/GWAS/')
    gwas_file = "mortality_COGs.txt"  # Replace with your actual file path
    fasta_file = "../Pangenome_results_HF/pan_genome_reference.fa"  # Replace with your actual file path
    output_file = "filtered_sequences_mortality.fasta"  # Output file for filtered sequences

    # Read the GWAS results and extract variant IDs
    variant_ids = read_table(gwas_file)

    # Extract sequences from the FASTA file based on the variant IDs
    sequences = extract_sequences(fasta_file, variant_ids)

    # Write the extracted sequences to a new FASTA file
    write_fasta(sequences, output_file)
    print(f"Filtered sequences written to {output_file}")
if __name__ == "__main__":
    main()
