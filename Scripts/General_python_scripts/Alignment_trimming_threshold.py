"""
Trim a multiple sequence alignment by removing columns that fall below a
minimum non-gap occupancy threshold.

For each column, the fraction of sequences with a non-gap character is
computed.  Columns where that fraction is below the specified threshold are
discarded.  The retained alignment is written to the output FASTA file and a
CSV recording removed column positions is saved alongside it.

Usage:
    python Alignment_trimming_threshold.py \
        --sequences <input.fasta> \
        --threshold <float, e.g. 0.8> \
        --output <output.fasta>

Arguments:
    --sequences  Path to the input FASTA alignment.
    --threshold  Minimum fraction (0–1) of non-gap characters required to
                 retain a column (e.g. 0.8 keeps columns with ≥80% non-gap).
    --output     Path for the trimmed output FASTA alignment.

Outputs:
    - <output>                       Trimmed alignment (FASTA)
    - removed_columns_threshold.csv  CSV of discarded column positions
"""

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd
import os

import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--sequences', required=True, type=str)
    parser.add_argument('--threshold', required=True, type=float)
    parser.add_argument('--output', required=True, type=str)
    params = parser.parse_args()


    # Load sequences
    sequences = params.sequences
    alignment = AlignIO.read(sequences, format="fasta")

    # Set the minimum percentage of non-gap characters to keep a column
    non_gap_threshold = params.threshold  # Adjust this value as needed

    # Calculate the minimum number of sequences without gaps required for each column
    min_non_gap_count = int(len(alignment) * non_gap_threshold)

    # Identify columns that meet the non-gap threshold
    alignment_length = alignment.get_alignment_length()
    columns_to_keep = []
    removed_columns_data = []

    for i in range(alignment_length):
        column_nucleotides = [record.seq[i] for record in alignment]
        non_gap_count = sum(1 for nuc in column_nucleotides if nuc != '-')

        if non_gap_count >= min_non_gap_count:
            columns_to_keep.append(i)
        else:
            removed_columns_data.append({
                'Column_Position': i,
                'Nucleotides': ''.join(column_nucleotides),
                'Non_Gap_Count': non_gap_count
            })

    # Convert removed columns to a Pandas DataFrame
    removed_columns_df = pd.DataFrame(removed_columns_data)

    # Extract and filter sequences based on retained columns
    filtered_sequences = []
    for record in tqdm(alignment):
        filtered_seq = ''.join(record.seq[i] for i in columns_to_keep)
        filtered_sequences.append(SeqRecord(Seq(filtered_seq), id=record.id, description=''))

    # Write the filtered sequences to output FASTA file
    output_path = params.output
    with open(output_path, 'w') as f:
        SeqIO.write(filtered_sequences, f, 'fasta')

    # Write the removed columns DataFrame to a CSV file
    csv_path = os.path.join(os.path.dirname(output_path), 'removed_columns_80thres.csv')
    removed_columns_df.to_csv(csv_path, index=False)

    print(f"Filtered sequences written to {output_path}")
    print(f"Removed columns log written to {csv_path}")
