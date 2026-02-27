"""
Build a multi-sample FASTA alignment from per-sample VCF-consensus FASTA files.

Iterates over all *.bam.fasta files in the specified directory, reads the
single consensus sequence from each, renames it using the filename stem, and
concatenates all records into one multi-sequence FASTA alignment file.

Usage:
    python Tree_VCF.py --directory <path/to/fasta_dir> --output <alignment.fasta>

Arguments:
    --directory  Directory containing *.bam.fasta consensus files.
    --output     Path for the combined output FASTA alignment.

Outputs:
    - <output>  Multi-sample FASTA alignment (one sequence per input file)
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse
from tqdm import tqdm

if '__main__' == __name__:

    parser = argparse.ArgumentParser()

    parser.add_argument('--directory', required=True, type=str)
    parser.add_argument('--output', required=True, type=str)
    params = parser.parse_args()

    os.chdir(params.directory)
    with open(params.output, 'a') as handle:
        for file in tqdm(os.listdir(params.directory)):
            if file.endswith('.bam.fasta'):
                fastafile = SeqIO.parse(file, 'fasta')
                for record in fastafile:
                    header = str(file).replace('.bam.fasta', '')
                    sequence = record.seq
                    seq_record = SeqRecord(sequence, id=header, description='')
                    SeqIO.write(seq_record, handle, 'fasta')
