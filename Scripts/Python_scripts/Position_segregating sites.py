from Bio import AlignIO
import numpy as np

# Define paths to alignment files
alignment_file = "/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Rlow_consensus_strains_onlyLucy.fasta"
lucy_alignment_file = "/home/albertotr/OneDrive/Data/Cambridge_Project/VCF_trees/Myco_py.fasta"

# Load alignments
alignment = AlignIO.read(alignment_file, "fasta")
alignment_lucy = AlignIO.read(lucy_alignment_file, "fasta")

def get_segregating_sites(alignment_input):
    """Identify segregating sites in an alignment."""
    segregating_sites = []
    for i in range(alignment_input.get_alignment_length()):
        column = set(record.seq[i] for record in alignment_input)
        if len(column) > 1:
            segregating_sites.append(i)
    return np.array(segregating_sites)

# Identify segregating sites
segregating_sites_1 = get_segregating_sites(alignment)
segregating_sites_2 = get_segregating_sites(alignment_lucy)

# Find common and unique elements
common_elements = np.intersect1d(segregating_sites_1, segregating_sites_2)
unique_to_list1 = np.setdiff1d(segregating_sites_1, segregating_sites_2)
unique_to_list2 = np.setdiff1d(segregating_sites_2, segregating_sites_1)

print("Common elements:", common_elements)
print("Unique to list 1:", unique_to_list1)
print("Unique to list 2:", unique_to_list2)
