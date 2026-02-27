"""
Split a multi-sequence FASTA alignment into individual per-sample FASTA files.

Reads a FASTA alignment and writes one file per sequence record to the specified
output directory.  A tab-delimited list file mapping sample IDs to output paths
is also created for downstream use (e.g., as Gubbins input list).

Usage:
    python Alignment2Fasta.py -i <input.fasta> -o <output_dir> -l <list_file.txt>

Arguments:
    -i / --input       Path to input FASTA alignment file.
    -o / --output_dir  Directory where individual FASTA files will be written.
    -l / --list_file   Path to output list file (TSV with ID and file path columns).

Outputs:
    - <output_dir>/<sample_id>.fasta  for each sequence in the input alignment
    - <list_file>  tab-separated list: sample_id <tab> absolute_path
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

from tqdm import tqdm
import argparse

# Load sequences
def main():
    parser = argparse.ArgumentParser(description="Divide a fasta file into multiple fasta files.")
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory for individual fasta files')
    parser.add_argument('-l', '--list_file', required=True, help='Output list file with paths to individual fasta files')
    args = parser.parse_args()

    sequences = args.input
    alignment = SeqIO.parse(sequences, format="fasta")

    # Change to the output directory
    os.chdir(args.output_dir)

    with open(args.list_file, 'a') as f:
        for record in tqdm(alignment):
            with open(str(record.id) + '.fasta', 'a') as handle:
                header = record.id
                sequence = record.seq
                seq_record = SeqRecord(sequence, id=header, description='')
                SeqIO.write(seq_record, handle, 'fasta')
                f.write(f"{header}\t{os.path.join(args.output_dir, header + '.fasta')}\n")

if __name__ == "__main__":
    main()
