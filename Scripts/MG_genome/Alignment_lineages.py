"""
Split a whole-dataset alignment into per-lineage sub-alignments.

Reads the Gubbins-masked, threshold-trimmed FASTA alignment for all samples
and partitions it into two sub-alignments based on pre-defined lineage
membership lists.  Lineage 2 samples are explicitly listed; all remaining
samples are assigned to Lineage 1.  Both sub-alignments are written as
separate FASTA files for downstream phylogenetic or comparative analyses.

Usage:
    Update the 'sequences', 'output_path1', and 'output_path2' variables at
    the top/bottom of the script, then run:
        python Alignment_lineages.py

Outputs:
    - VA94_lineage1_60threshold.masked.fasta  Lineage 1 sub-alignment
    - VA94_lineage2_60threshold.masked.fasta  Lineage 2 sub-alignment
"""

from Bio import AlignIO

# Load sequences
sequences ='/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_consensus_all_trimmed_60threshold.fasta'
alignment = AlignIO.read(sequences, format="fasta")

# Lineage 2
lineage2 = [
'A090809_2009',
'G_2015',
'L_2015',
'J_2015',
'A01546_2007',
'A01505_2007',
'A01510_2007',
'A565_2002',
'A2261_2000',
'A506_2002',
'A454_2001',
'A442_2001',
'A471_2001',
'A267_2000',
'A3611_2001',
'A509_2002',
'A277_2001',
'BP421_2002',
'A933_2001',
'A297_2001',
'A371_2001',
'A877_2003',
'MG31_AL_11_2011',
'A020_2011',
'OY79_2013',
'WU47_2013',
'A_2015',
'Q_2015',
'Black_2012',
'GW77_2013',
'E054_2013',
'KB64_2013',
'A1_2014',
'MG27_AL_11_2011',
'A046_2011',
'A3278_2012',
'A3240_2012',
'A3244_2012',
'A3225_2012',
'A2486_2012',
'A304_2003',
'A195_2003',
'A012_2011',
'MG28_AL_11_2011',
'A044_2011',
'A001_2011',
'MG22_AL_11_2011',
'A004_2011',
'A3316_2012',
'Blue_2012',
'A021_2011',
'MG23_AL_11_2011',
'A006_2011',
'MG25_AL_11_2011',
'MG32_AL_11_2011',
'AMGO_2011',
'MG29_AL_11_2011',
'A013_2011',
'MG26_AL_11_2011',
'A018_2011',
'MG30_AL_11_2011',
'NC08_2008.031-4-3P',
'NC06_2006.080-5-2P'
            ]

# Separate lineage 2 sequences into a new alignment
lineage1_sequences = [record for record in alignment if record.id not in lineage2]
lineage1_alignment = AlignIO.MultipleSeqAlignment(lineage1_sequences)
lineage2_sequences = [record for record in alignment if record.id in lineage2]
lineage2_alignment = AlignIO.MultipleSeqAlignment(lineage2_sequences)

#%% Write both alignments to output FASTA files
output_path1 = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_lineage1_60threshold.masked.fasta'
output_path2 = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_lineage2_60threshold.masked.fasta'

with open(output_path1, 'w') as f:
    AlignIO.write(lineage1_alignment, f, 'fasta')
with open(output_path2, 'w') as f:
    AlignIO.write(lineage2_alignment, f, 'fasta')