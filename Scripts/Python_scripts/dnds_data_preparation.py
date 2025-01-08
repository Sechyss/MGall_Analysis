from BCBio import GFF
from Bio import AlignIO, Align
from Bio.Seq import Seq

# Initialize dictionaries to store gene data
collector_dict = {}
for record in GFF.parse('/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCA_000286675.1/genomic.gff'):
    for genes in record.features:
        if genes.id == '':
            continue
        else:
            id_index = str(genes.id).replace('gene-', '')
            start = int(genes.location.start)
            end = int(genes.location.end)
            collector_dict.update({id_index: [start, end]})

# Remove irrelevant entries
collector_dict_filtered = {key: value for key, value in collector_dict.items() if 'CP003506.1' not in key}

# Read the alignment
alignment = AlignIO.read('/home/albertotr/OneDrive/Data/Cambridge_project/Mapped_output_VA94_7994_1_7P/VA94_consensus_spns.masked.aln', 'fasta')

# Flatten gene ranges into a single list of positions to keep
positions_to_keep = set()
for start, end in collector_dict_filtered.values():
    positions_to_keep.update(range(start, end + 1))

# Identify columns with no gaps in any sequence
alignment_length = alignment.get_alignment_length()
columns_to_keep = []
for i in range(alignment_length):
    column_nucleotides = [record.seq[i] for record in alignment]
    if any(nuc != '-' for nuc in column_nucleotides):
        columns_to_keep.append(i)

# Create a new alignment with only the columns within the gene regions and no gap-only columns
trimmed_records = []
for record in alignment:
    trimmed_seq = ''.join([record.seq[i] for i in columns_to_keep if i + 1 in positions_to_keep])
    trimmed_records.append(record[:0])  # Copy the record structure
    trimmed_records[-1].seq = Seq(trimmed_seq)

# Create a new alignment object
trimmed_alignment = Align.MultipleSeqAlignment(trimmed_records)

# Save the trimmed alignment
output_path = '/home/albertotr/OneDrive/Data/Cambridge_project/Mapped_output_VA94_7994_1_7P/VA94_consensus_spns.masked_trimmed.aln'
AlignIO.write(trimmed_alignment, output_path, 'fasta')