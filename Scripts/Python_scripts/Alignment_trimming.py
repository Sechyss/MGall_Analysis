from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Load sequences
sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Rlow_consensus_spns.masked.aln'
alignment = AlignIO.read(sequences, format="fasta")

# Identify columns with no gaps in any sequence
alignment_length = alignment.get_alignment_length()
columns_to_keep = [
    i for i in range(alignment_length)
    if all(record.seq[i] != '-' for record in alignment)
]

# Extract sequences based on gap-free columns
filtered_sequences = []
for record in tqdm(alignment):
    filtered_seq = ''.join(record.seq[i] for i in columns_to_keep)
    filtered_sequences.append(SeqRecord(Seq(filtered_seq), id=record.id, description=''))

# Write the filtered sequences to output FASTA file
output_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Rlow_consensus_snps_trimmed_nogaps.masked.fasta'
with open(output_path, 'w') as f:
    SeqIO.write(filtered_sequences, f, 'fasta')
