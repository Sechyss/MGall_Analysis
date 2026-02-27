"""
Extract genome FASTA sequences from Prokka GFF files.

Reads the embedded FASTA section from each GFF file in an input directory,
renames the sequence using a pre-built name-replacement dictionary (Lucy
sample IDs → standardised identifiers), and writes individual .fna files
to an output directory.

Usage:
    Update 'input_folder', 'output_folder', and the path to the replacement
    pickle file at the bottom of this script, then run:
        python GFFs2FNA.py

Outputs:
    - <output_folder>/<new_id>.fna  one per GFF with a matching replacement key
"""

import os
import pickle
from Bio import SeqIO

# Dictionary to translate the name of the file
# Load the dictionary with the replacements
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'rb') as handle:
    sorted_replacements = pickle.load(handle)

def extract_fna_from_gff(gff_file, output_folder):
    file_base_name = os.path.basename(gff_file).replace(".gff", "")
    if file_base_name in sorted_replacements:
        new_id = sorted_replacements[file_base_name]
        with open(gff_file, 'r') as file:
            genome_sequence = ""
            for line in file:
                if line.startswith("##FASTA"):
                    break
            for record in SeqIO.parse(file, "fasta"):
                genome_sequence += str(record.seq)
            output_file = os.path.join(output_folder, f"{new_id}.fna")
            with open(output_file, 'w') as out_file:
                out_file.write(f">{new_id}\n{genome_sequence}\n")

def process_gff_folder(input_folder, output_folder):
    for filename in os.listdir(input_folder):
        if filename.endswith(".gff"):
            gff_file = os.path.join(input_folder, filename)
            extract_fna_from_gff(gff_file, output_folder)

if __name__ == "__main__":
    input_folder = "/mnt/c/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/GFFs_genomes/HF_gffs/"
    output_folder = "/mnt/c/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Pangenome_bins/"
    process_gff_folder(input_folder, output_folder)
