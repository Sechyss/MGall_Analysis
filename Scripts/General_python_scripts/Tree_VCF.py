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
