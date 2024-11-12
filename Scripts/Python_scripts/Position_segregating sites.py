import os
from Bio import AlignIO
import pandas as pd
import numpy as np
import pysam

# Define paths to alignment files
alignment_file = "/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Rlow_consensus_strains_onlyLucy.fasta"
lucy_alignment_file = "/home/albertotr/OneDrive/Data/Cambridge_Project/VCF_trees/Myco_py_noref.fasta"

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
    nucleotides_1 = set(record.seq[pos].upper() for record in alignment if record.seq[pos] != '-')
    nucleotides_2 = set(record.seq[pos].upper() for record in alignment_lucy if record.seq[pos] != '-')

    # Skip if either set is empty or if both sets have the same nucleotides
    if not nucleotides_1 or not nucleotides_2 or nucleotides_1 == nucleotides_2:
        continue

    # Skip if there is any common letter between the two sets
    if nucleotides_1 & nucleotides_2:
        continue

    # Add to list if nucleotides are different
    diff_data.append({
        "Position": pos,
        "Alignment 1 Nucleotides": ''.join(nucleotides_1),
        "Alignment 2 Nucleotides": ''.join(nucleotides_2)
    })

diff_df = pd.DataFrame(diff_data)

del segregating_sites_1, segregating_sites_2
# Complement the table with the BAM Files
# List of positions to check in BAM files
positions = diff_df["Position"].tolist()

# Directory containing BAM files
bam_dir = "/home/albertotr/downloads/BAM/BAM/"

#%%  Complement with BAM files

# Loop through each BAM file and check sequences at differing positions
for file in os.listdir(bam_dir):
    if file.endswith(".bam"):
        bamfile_path = os.path.join(bam_dir, file)
        bamfile = pysam.AlignmentFile(bamfile_path, "rb")

        # Assuming single chromosome reference; adjust if necessary
        chromosome = bamfile.references[0]

        # Store sequences for each position
        sequences = []
        for pos_nucleotide in positions:
            bases = []
            for pileupcolumn in bamfile.pileup(chromosome, pos_nucleotide, pos_nucleotide + 1, truncate=True):
                if pileupcolumn.pos == pos_nucleotide:  # BAM is 0-based, so this matches our positions
                    bases = [
                        pileupread.alignment.query_sequence[pileupread.query_position]
                        for pileupread in pileupcolumn.pileups
                        if pileupread.query_position is not None
                    ]
                    break
            sequences.append(''.join(bases))

        # Add sequences as a new column in diff_df
        bamfile_label = os.path.splitext(file)[0]  # Name column based on BAM filename
        diff_df[bamfile_label] = sequences

        bamfile.close()
diff_df.to_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Test.csv', index=False)

del pileupcolumn

#%% Alternative code to check the position of the reads

df_pos = pd.DataFrame(columns=['Read_name', 'BAM_File', 'Position', 'Location'])

for file in os.listdir(bam_dir):
    if file.endswith(".bam"):
        bamfile_path = os.path.join(bam_dir, file)
        bamfile = pysam.AlignmentFile(bamfile_path, "rb")

        # Assuming single chromosome reference; adjust if necessary
        chromosome = bamfile.references[0]
        for pos_nucleotide in positions:
            for pileupcolumn in bamfile.pileup(chromosome, pos_nucleotide, pos_nucleotide +1, truncate=True):  # BAM is 0-based
                if pileupcolumn.pos == pos_nucleotide:  # Match 0-based BAM with your 1-based position
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.query_position is None:
                            continue

                        read_position = pileupread.query_position
                        read_length = pileupread.alignment.query_length

                        # Check if nucleotide is at the beginning, end, or middle of the read
                        if 0 <= read_position < 10:
                            location = "beginning"
                        elif read_length - 10 <= read_position < read_length:
                            location = "end"
                        else:
                            location = "middle"

                        row_to_add = pd.DataFrame([{'Read_name':pileupread.alignment.query_name, 'BAM_File': os.path.splitext(file)[0],
                                      'Position': pos_nucleotide, 'Location': location}])
                        # Print information
                        df_pos = pd.concat([df_pos, row_to_add], ignore_index=True)
        bamfile.close()
df_pos.to_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/Test_pos.csv', index=False)
del pileupread, pileupcolumn