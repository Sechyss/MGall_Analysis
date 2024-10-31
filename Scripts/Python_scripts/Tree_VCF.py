from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

from tqdm import tqdm

os.chdir('/home/albertotr/OneDrive/Data/'
         'Cambridge_Project/Mapped_output_Rlow/Only_SNPs/')
with open('Rlow_consensus_strains_onlyLucy.fasta', 'a') as handle:
    for file in tqdm(os.listdir('/home/albertotr/OneDrive/Data/'
                                'Cambridge_Project/Mapped_output_Rlow/Only_SNPs/')):
        if file.endswith('.bam.fasta'):
            fastafile = SeqIO.parse(file, 'fasta')
            for record in fastafile:
                header = str(file).replace('.bam.fasta', '')
                sequence = record.seq
                seq_record = SeqRecord(sequence, id=header, description='')
                SeqIO.write(seq_record, handle, 'fasta')
