"""
Filter a multiple sequence alignment to parsimony-informative sites only.

Identifies columns that are parsimony-informative (at least two distinct
character states, each present in ≥2 sequences, ignoring gaps) and writes
a new alignment containing only those columns.

Usage:
    Update the 'sequences' and 'output_alignment' path variables at the bottom
    of the script, then run:
        python Alignment_trimming_informativesites.py

Outputs:
    - Filtered FASTA alignment retaining only parsimony-informative columns
"""

from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from collections import Counter

from Bio.Seq import Seq


def find_parsimony_informative_sites(alignment_file, format_file="fasta"):
    """
    Identify parsimony-informative sites in an alignment.

    Parameters:
    - alignment_file: Path to the alignment file.
    - format: Format of the alignment file (e.g., 'fasta', 'phylip').

    Returns:
    - List of column indices that are parsimony-informative.
    """
    # Load the alignment
    alignment = AlignIO.read(alignment_file, format_file)
    informative_sites = []

    for col_idx in range(alignment.get_alignment_length()):
        column = [record.seq[col_idx] for record in alignment]
        counts = Counter(column)

        # Ignore gaps
        counts.pop('-', None)

        # Check for parsimony-informativeness
        if len(counts) >= 2 and sum(1 for v in counts.values() if v >= 2) >= 2:
            informative_sites.append(col_idx)

    return informative_sites


sequences = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_all_gubbins_run.masked.aln'
informative_sites_bp = find_parsimony_informative_sites(sequences)


def calculate_column_variability(alignment_file, format_file="fasta"):
    """
    Calculate variability for each column in an alignment.

    Parameters:
    - alignment_file: Path to the alignment file.
    - format: Format of the alignment file (e.g., 'fasta', 'phylip').

    Returns:
    - List of variability scores (proportion of unique characters per column).
    """
    alignment = AlignIO.read(alignment_file, format_file)
    variability = []

    for col_idx in range(alignment.get_alignment_length()):
        column = [record.seq[col_idx] for record in alignment]
        unique_chars = set(column) - {'-'}  # Exclude gaps
        variability.append(len(unique_chars) / len(column))

    return variability


# Example usage
variability_scores = calculate_column_variability(sequences)

# Filter columns based on variability

def filter_columns_by_informative_sites(alignment_file, output_file, informative_sites, format_file="fasta"):
    """
    Retain only parsimony-informative columns in the alignment.

    Parameters:
    - alignment_file: Path to the input alignment file.
    - output_file: Path to save the filtered alignment.
    - informative_sites: List of column indices to retain (0-based).
    - format_file: Format of the alignment file (default: 'fasta').
    """
    # Load the alignment
    alignment = AlignIO.read(alignment_file, format_file)
    filtered_records = []

    for record in alignment:
        # Filter sequence based on informative columns
        filtered_seq = ''.join(record.seq[i] for i in informative_sites)
        record.seq = Seq(filtered_seq)
        filtered_records.append(record)

    # Convert filtered records to MultipleSeqAlignment
    filtered_alignment = MultipleSeqAlignment(filtered_records)

    # Write the filtered alignment to a file
    AlignIO.write(filtered_alignment, output_file, format_file)


# Example usage
output_alignment = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VA94_consensus_all_trimmed_informative_sites.fasta'
filter_columns_by_informative_sites(sequences, output_alignment, informative_sites_bp)
