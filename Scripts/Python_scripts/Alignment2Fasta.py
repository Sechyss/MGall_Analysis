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
