"""
Verify and visualise gene context (synteny neighbourhood) around target loci.

Parses GFF annotation files and a presence/absence summary table to build
neighbourhood diagrams for a set of query genes.  For each query, the flanking
coding sequences are extracted, mapped through a replacement-value dictionary
to standardised product/cluster names, and drawn as a linear genomic feature
diagram using dna_features_viewer.

Usage:
    Update the path variables and target gene list inside the script, then run:
        python GeneContext_verification.py

Outputs:
    - Per-gene PNG/SVG context diagrams (one per query locus)
"""

import os

import numpy as np
from tqdm import tqdm

import pandas as pd
from Bio import SeqIO
from BCBio import GFF
import matplotlib.patches as mpatches
from dna_features_viewer import GraphicFeature, GraphicRecord


def replace_value(cell, dictionary_info):
    # Split the cell by ';' to handle multiple values
    values_to_add = cell.split(';')
    new_values = []
    for value in values_to_add:
        replaced = False
        for key, value_dict in dictionary_info.items():
            if value in value_dict:
                new_values.append(str(value_dict[value]))
                replaced = True
                break
        if not replaced:
            new_values.append(value)
    return ';'.join(new_values)


def flatten_list(nested_list):
    flat_list = []
    for sublist in nested_list:
        for item in sublist:
            flat_list.append(item)
    return flat_list


fastafile = SeqIO.parse('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                        'pangenome_results_filtered/Genes_to_study_post2007_HF.faa', 'fasta')
list_fasta_header = []
for record in fastafile:
    id_to_add = record.id
    list_fasta_header.append(id_to_add)

Presence_absence = pd.read_csv('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                               'pangenome_results_filtered/gene_presence_absence.csv')

potential_phage = flatten_list(pd.read_excel('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                                             'pangenome_results_filtered/Genestostudy_HF_post2007_hitdata.xlsx',
                                             sheet_name='Potential_phage_genes').values.tolist())

filtered_table = Presence_absence[Presence_absence['Gene'].isin(potential_phage)].drop(
    columns=['Non-unique Gene name', 'Annotation']).set_index(['Gene'])
filtered_table.dropna(axis=1, how='all', inplace=True)

collector_dict = {}
strand_dict = {}
for file in tqdm(os.listdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                            'GFFs_genomes/filtered_genomes/')):
    filename = file.split('.')[0]
    genes_dict = {}
    gene_strand = {}
    for record in GFF.parse('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                            'GFFs_genomes/filtered_genomes/' + file):
        for genes in record.features:
            if genes.id == '':
                continue
            else:
                id_index = genes.id
                start = int(genes.location.start)
                end = int(genes.location.end)
                genes_dict.update({id_index: [start, end]})
                gene_strand.update({id_index: genes.location.strand})
    collector_dict.update({filename: genes_dict})
    strand_dict.update({filename: gene_strand})

new_table = filtered_table.map(lambda cell: replace_value(str(cell), collector_dict))
new_table = new_table.map(lambda cell: np.nan if '[' not in str(cell) else cell)

new_table_strand = filtered_table.map(lambda cell: replace_value(str(cell), strand_dict))
new_table_strand = new_table_strand.map(lambda cell: np.nan if cell != int else cell)

for col in new_table:
    series = new_table[col].dropna()
    ends = []
    coordinate_features = []
    for index, coordinates in series.items():
        if coordinates != 'nan':
            if ';' in coordinates:
                first, second = coordinates.replace('[', '').replace(']', '').split(';')
                first_start = int(first.split(',')[0])
                first_end = int(first.split(',')[1])
                second_start = int(second.split(',')[0])
                second_end = int(second.split(',')[1])
                features_to_add_1 = GraphicFeature(start=first_start, end=first_end, color="#ffd700",
                                                   label=index)
                features_to_add_2 = GraphicFeature(start=second_start, end=second_end,
                                                   color="#ffd700", label=index)
                coordinate_features.append(features_to_add_1)
                coordinate_features.append(features_to_add_2)
                ends.append(first_end)
                ends.append(second_end)

            else:
                pos = coordinates.replace('[', '').replace(']', '').split(',')
                first_start = int(pos[0])
                first_end = int(pos[1])
                features_to_add = GraphicFeature(start=first_start, end=first_end, color="#ffd700",
                                                 label=index)
                coordinate_features.append(features_to_add)
                ends.append(first_end)

    max_end = max(ends)
    record = GraphicRecord(sequence_length=max_end, features=coordinate_features)
    ax, _ = record.plot()
    patch0 = mpatches.Patch(color='#ffd700', label='Potential phage genes')
    legend = ax.legend(handles=[patch0], loc='upper right', fontsize=10,
                       title=str(col), title_fontsize=15)
    legend.get_frame().set_alpha(None)
    ax.figure.savefig('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                      'pangenome_results_filtered/Phage_gene_location/' + str(col) + '.pdf',
                      bbox_inches='tight')
