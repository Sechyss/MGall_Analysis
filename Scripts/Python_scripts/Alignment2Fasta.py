from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

from tqdm import tqdm

# Load sequences
sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_consensus_all_60thres.fasta'
alignment = SeqIO.parse(sequences, format="fasta")

#%% Divide into multiple fasta files

path_2_file = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/Trimmed_fasta/'
os.chdir(path_2_file)

with open('/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/List_files_all_VA94-Trimmed.txt', 'a') as f:
    for record in tqdm(alignment):
        with open(str(record.id)+'.fasta', 'a') as handle:
            header = record.id
            sequence = record.seq
            seq_record = SeqRecord(sequence, id=header, description='')
            SeqIO.write(seq_record, handle, 'fasta')
            f.write(f"{header}\t{path_2_file + header+'.fasta'}\n")
