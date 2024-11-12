from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Load sequences
sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Rlow_consensus_spns.masked.aln'
alignment = AlignIO.read(sequences, format="fasta")

# Set the minimum percentage of non-gap characters to keep a column
non_gap_threshold = 0.8  # Adjust this value as needed

# Calculate the minimum number of sequences without gaps required for each column
min_non_gap_count = int(len(alignment) * non_gap_threshold)

# Identify columns with non-gap counts meeting the threshold
alignment_length = alignment.get_alignment_length()
columns_to_keep = [
    i for i in range(alignment_length)
    if sum(1 for record in alignment if record.seq[i] != '-') >= min_non_gap_count
]

# Extract and filter sequences based on retained columns
filtered_sequences = []
for record in tqdm(alignment):
    filtered_seq = ''.join(record.seq[i] for i in columns_to_keep)
    filtered_sequences.append(SeqRecord(Seq(filtered_seq), id=record.id, description=''))

# Write the filtered sequences to output FASTA file
output_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Rlow_consensus_onlyLucy_trimmed_threshold80.fasta'
with open(output_path, 'w') as f:
    SeqIO.write(filtered_sequences, f, 'fasta')
