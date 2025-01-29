from BCBio import GFF
from Bio import AlignIO, Align
from Bio.Seq import Seq

# Initialize dictionaries to store gene data
collector_dict = {}
for record in GFF.parse('/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCA_000286675.1/genomic.gff'):
    for genes in record.features:
        if genes.id == '' or 'CP003506.1' in genes.id:
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

# Identify columns with no gaps in any sequence and are part of a gene
alignment_length = alignment.get_alignment_length()
columns_to_keep = []
for i in range(alignment_length):
    if i + 1 in positions_to_keep:  # Check if the position is part of a gene
        column_nucleotides = [record.seq[i] for record in alignment]
        if any(nuc != '-' for nuc in column_nucleotides):
            columns_to_keep.append(i)

# Create a new alignment with only the columns within the gene regions and no gap-only columns
trimmed_records = []
for record in alignment:
    trimmed_seq = ''.join([record.seq[i] for i in columns_to_keep])
    trimmed_records.append(record[:0])  # Copy the record structure
    trimmed_records[-1].seq = Seq(trimmed_seq)

# Create a new alignment object
trimmed_alignment = Align.MultipleSeqAlignment(trimmed_records)

# Save the trimmed alignment
output_path = '/home/albertotr/OneDrive/Data/Cambridge_project/Mapped_output_VA94_7994_1_7P/dnds/VA94_test.masked.aln'
AlignIO.write(trimmed_alignment, output_path, 'fasta')

#%% Create individual alignments for each gene
# Loop through each gene in the filtered collector dictionary
for gene_id, (start, end) in collector_dict_filtered.items():
    # Create a range of positions for the current gene
    gene_positions = range(start, end + 1)
    # Identify columns in the alignment that correspond to the gene positions
    gene_columns = [i for i in columns_to_keep if i + 1 in gene_positions]
    # Initialize a list to store the gene-specific records
    gene_records = []
    # Loop through each record in the alignment
    for record in alignment:
        # Extract the sequence for the current gene based on the identified columns
        gene_seq = ''.join([record.seq[i] for i in gene_columns])
        # Copy the record structure and assign the extracted gene sequence
        gene_records.append(record[:0])
        gene_records[-1].seq = Seq(gene_seq)
    # Create a new alignment object for the current gene
    gene_alignment = Align.MultipleSeqAlignment(gene_records)
    # Define the output path for the gene-specific alignment file
    gene_output_path = f'/home/albertotr/OneDrive/Data/Cambridge_project/Mapped_output_VA94_7994_1_7P/dnds/genes/{gene_id}.aln'
    # Save the gene-specific alignment
    AlignIO.write(gene_alignment, gene_output_path, 'fasta')

    # Translate the gene sequences to protein sequences using genetic table 4
    protein_records = []
    for record in gene_records:
        protein_seq = record.seq.translate(table=4)
        protein_record = record.__class__(record.id, protein_seq, record.description)  # Copy the record structure
        protein_record.seq = protein_seq
        protein_records.append(protein_record)
    
    # Create a new alignment object for the translated proteins
    protein_alignment = Align.MultipleSeqAlignment(protein_records)
    # Define the output path for the protein-specific alignment file
    protein_output_path = f'/home/albertotr/OneDrive/Data/Cambridge_project/Mapped_output_VA94_7994_1_7P/dnds/proteins/{gene_id}.aln'
    # Save the protein-specific alignment
    AlignIO.write(protein_alignment, protein_output_path, 'fasta')