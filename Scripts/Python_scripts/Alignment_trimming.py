import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Load sequences
sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Edited_Rlow_consensus_strains_onlyLucy.fasta'
alignment = AlignIO.read(sequences, format="fasta")

# Create column coordinates based on the length of the sequence
coordinates = list(range(1, len(alignment[0].seq) + 1))

# Build DataFrame directly
data = {record.id: list(record.seq) for record in tqdm(alignment)}
df = pd.DataFrame(data, index=coordinates).T

# Remove columns with more than 60% gaps
threshold = len(df) * 0.6
df2 = df.loc[:, (df == '-').sum(axis=0) <= threshold]

# Combine each row into a single string for FASTA output
df2['Sequence'] = df2.apply(lambda row: ''.join(row.astype(str)), axis=1)

# Write the trimmed sequences to FASTA file in one go
output_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Edited_Rlow_consensus_strains_onlyLucy_trimmed.fasta'
seq_records = [SeqRecord(seq=Seq(row['Sequence']), id=str(index), description='') for index, row in df2.iterrows()]
with open(output_path, 'w') as f:
    SeqIO.write(seq_records, f, 'fasta')
