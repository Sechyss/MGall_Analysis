import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Edited_Rlow_consensus_strains_onlyLucy.fasta'

alignment = AlignIO.read(sequences, format="fasta")
counter = 1
coordinates = []
for record in alignment:
        sequence = record.seq
        for aminoacid in sequence:
            coordinates.append(counter)
            counter += 1
        break

df = pd.DataFrame(index=[records.id for records in alignment], columns=coordinates)

for i, col in tqdm(enumerate(alignment)):
    aligment_id = col.id
    all_seq = col.seq
    for index, aminoacid in enumerate(all_seq):
        df.loc[aligment_id, coordinates[index]] = aminoacid

# Calculate the threshold for 60% of rows
threshold = len(df) * 0.6

# Drop columns where '-' appears in more than 60% of the cells
df2 = df.drop(columns=[col for col in df.columns if (df[col] == '-').sum() > threshold])

# Combine all data in each row into a single string
df2['Sequence'] = df2.apply(lambda row: ' '.join(row.astype(str)), axis=1)

with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Edited_Rlow_consensus_strains_onlyLucy_trimmed.fasta') as f:
    for index, row in df2.iterrows():
        sequence = row['Sequence']
        sequence_id =  index
        seq_record = SeqRecord(seq=Seq(sequence), id=str(sequence_id), description='')
        SeqIO.write(seq_record, f, 'fasta')
