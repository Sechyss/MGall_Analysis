"""
Convert gene-level .aln.fas alignment files to tree-ready FASTA format.

Iterates over all *.aln.fas files in the specified directory, reformats the
sequence headers to '<sample_id>_<gene_name>' (splitting on ';' and removing
problematic '~' characters), and writes the renamed sequences to new .fasta
files in the same directory.

Usage:
    python Preparation_trees.py --directory <path/to/alignment_dir>

Arguments:
    --directory  Directory containing *.aln.fas alignment files.

Outputs:
    - <gene>.fasta  Renamed alignment file for each *.aln.fas input
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--directory', required=True, type=str)
    params = parser.parse_args()

    os.chdir(params.directory)
    for file in os.listdir(params.directory):
        if not file.endswith('.aln.fas'):
            continue
        else:
            fastafile = SeqIO.parse(file, 'fasta')
            with open(str(file).replace('.aln.fas', '.fasta').replace('~', '-'), 'a') as handle:
                for record in fastafile:
                    header = str(record.id).split(';')[0]+'_'+str(file).replace('.aln.fas', '').replace('~', '-')
                    sequence = record.seq
                    seq_record = SeqRecord(sequence, id=header, description='')
                    SeqIO.write(seq_record, handle, 'fasta')
