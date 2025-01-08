from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd

# Load sequences
sequences ='/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_all_gubbins_run.masked.aln'
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
output_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_consensus_all_trimmed_nogaps.masked.fasta'
with open(output_path, 'w') as f:
    SeqIO.write(filtered_sequences, f, 'fasta')

# Write the removed columns DataFrame to a CSV file
csv_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/removed_columns.csv'
removed_columns_df.to_csv(csv_path, index=False)
