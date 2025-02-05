from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Load sequences
sequences ='/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_all_gubbins_run.masked.aln'
alignment = AlignIO.read(sequences, format="fasta")

# Lineage 2
lineage2 = ['A006_2011',
            'A013_2011',
            'A018_2011',
            'A021_2011',
            'AMGO_2011',
            'MG23_AL_11_2011',
            'MG25_AL_11_2011',
            'MG26_AL_11_2011',
            'MG29_AL_11_2011',
            'MG30_AL_11_2011',
            'MG32_AL_11_2011'
            ]

# Separate lineage 2 sequences into a new alignment
lineage1_sequences = [record for record in alignment if record.id not in lineage2]
lineage1_alignment = AlignIO.MultipleSeqAlignment(lineage1_sequences)
lineage2_sequences = [record for record in alignment if record.id in lineage2]
lineage2_alignment = AlignIO.MultipleSeqAlignment(lineage2_sequences)

# Write both alignments to output FASTA files
output_path1 = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_lineage1.masked.fasta'
output_path2 = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_lineage2.masked.fasta'

with open(output_path1, 'w') as f:
    AlignIO.write(lineage1_alignment, f, 'fasta')
with open(output_path2, 'w') as f:
    AlignIO.write(lineage2_alignment, f, 'fasta')