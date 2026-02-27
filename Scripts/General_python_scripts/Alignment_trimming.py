"""
Remove gap-only (all-gap) columns from a multiple sequence alignment.

Reads a FASTA alignment and discards any column position where every sequence
carries a gap character ('-').  The filtered alignment is written to a new FASTA
file and a CSV recording the removed column positions and nucleotide composition
is saved alongside the output.

Usage:
    python Alignment_trimming.py --sequences <input.fasta> --output <output.fasta>

Arguments:
    --sequences  Path to the input FASTA alignment.
    --output     Path for the trimmed output FASTA alignment.

Outputs:
    - <output>                   Trimmed alignment (FASTA)
    - removed_columns_nogaps.csv CSV of discarded column positions and nucleotides
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
    parser.add_argument('--output', required=True, type=str)
    params = parser.parse_args()

    # Load sequences
    sequences =params.sequences
    alignment = AlignIO.read(sequences, format="fasta")

    # Identify columns with no gaps in any sequence
    alignment_length = alignment.get_alignment_length()
    columns_to_keep = []
    removed_columns_data = []

    for i in range(alignment_length):
        column_nucleotides = [record.seq[i] for record in alignment]
        if all(nuc != '-' for nuc in column_nucleotides):
            columns_to_keep.append(i)
        else:
            removed_columns_data.append({
                'Column_Position': i,
                'Nucleotides': ''.join(column_nucleotides)
            })

    # Convert removed columns to a Pandas DataFrame
    removed_columns_df = pd.DataFrame(removed_columns_data)

    # Extract sequences based on gap-free columns
    filtered_sequences = []
    for record in tqdm(alignment):
        filtered_seq = ''.join(record.seq[i] for i in columns_to_keep)
        filtered_sequences.append(SeqRecord(Seq(filtered_seq), id=record.id, description=''))

    # Write the filtered sequences to output FASTA file
    output_path = params.output
    with open(output_path, 'w') as f:
        SeqIO.write(filtered_sequences, f, 'fasta')

    # Write the removed columns DataFrame to a CSV file
    csv_path = os.path.join(os.path.dirname(output_path), 'removed_columns_nogaps.csv')
    removed_columns_df.to_csv(csv_path, index=False)
