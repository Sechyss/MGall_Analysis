from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

from tqdm import tqdm

# Load sequences
sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/VA94_consensus_onlyLucy_trimmed_threshold60.fasta'
alignment = SeqIO.parse(sequences, format="fasta")

#%% Divide into multiple fasta files

path_2_file = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/Trimmed_fasta/'
os.chdir(path_2_file)

with open('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/List_files_VA94-Trimmed.txt', 'a') as f:
    for record in tqdm(alignment):
        with open(str(record.id)+'.fasta', 'a') as handle:
            header = record.id
            sequence = record.seq
            seq_record = SeqRecord(sequence, id=header, description='')
            SeqIO.write(seq_record, handle, 'fasta')
            f.write(f"{header}\t{path_2_file + header+'.fasta'}\n")
