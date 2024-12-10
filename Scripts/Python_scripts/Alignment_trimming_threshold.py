from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pandas as pd

# Load sequences
sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/VA94_consensus_spns.masked.aln'
alignment = AlignIO.read(sequences, format="fasta")

# Set the minimum percentage of non-gap characters to keep a column
non_gap_threshold = 0.6  # Adjust this value as needed

# Calculate the minimum number of sequences without gaps required for each column
min_non_gap_count = int(len(alignment) * non_gap_threshold)

# Identify columns that meet the non-gap threshold
alignment_length = alignment.get_alignment_length()
columns_to_keep = []
removed_columns_data = []

for i in range(alignment_length):
    column_nucleotides = [record.seq[i] for record in alignment]
    non_gap_count = sum(1 for nuc in column_nucleotides if nuc != '-')

    if non_gap_count >= min_non_gap_count:
        columns_to_keep.append(i)
    else:
        removed_columns_data.append({
            'Column_Position': i,
            'Nucleotides': ''.join(column_nucleotides),
            'Non_Gap_Count': non_gap_count
        })

# Convert removed columns to a Pandas DataFrame
removed_columns_df = pd.DataFrame(removed_columns_data)

# Extract and filter sequences based on retained columns
filtered_sequences = []
for record in tqdm(alignment):
    filtered_seq = ''.join(record.seq[i] for i in columns_to_keep)
    filtered_sequences.append(SeqRecord(Seq(filtered_seq), id=record.id, description=''))

# Write the filtered sequences to output FASTA file
output_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/VA94_consensus_onlyLucy_trimmed_threshold60.fasta'
with open(output_path, 'w') as f:
    SeqIO.write(filtered_sequences, f, 'fasta')

# Write the removed columns DataFrame to a CSV file
csv_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/removed_columns_threshold60.csv'
removed_columns_df.to_csv(csv_path, index=False)

print(f"Filtered sequences written to {output_path}")
print(f"Removed columns log written to {csv_path}")
