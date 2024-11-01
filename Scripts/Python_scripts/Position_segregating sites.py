from Bio import AlignIO
import pandas as pd
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

# Find differing elements (symmetric difference)
different_sites = np.setxor1d(segregating_sites_1, segregating_sites_2)

# Create a DataFrame for differing sites
diff_data = []
for pos in different_sites:
    nucleotides_1 = set(record.seq[pos] for record in alignment if record.seq[pos] != '-')
    nucleotides_2 = set(record.seq[pos] for record in alignment_lucy if record.seq[pos] != '-')

    # Skip if either set is empty or if both sets have the same nucleotides
    if not nucleotides_1 or not nucleotides_2 or nucleotides_1 == nucleotides_2:
        continue

    # Add to list if nucleotides are different
    diff_data.append({
        "Position": pos,
        "Alignment 1 Nucleotides": ''.join(nucleotides_1),
        "Alignment 2 Nucleotides": ''.join(nucleotides_2)
    })

diff_df = pd.DataFrame(diff_data)
