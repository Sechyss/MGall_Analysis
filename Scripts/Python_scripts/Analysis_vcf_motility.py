import os

import pandas as pd
import pysam
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Initialize dictionaries to store gene data
collector_dict = {}
gene_names = {}

# Parse the GFF file to populate dictionaries
for record in GFF.parse('/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCF_000286675.1/genomic.gff'):
    for genes in record.features:
        if genes.id == '':
            continue
        else:
            id_index = str(genes.id).replace('gene-', '')
            start = int(genes.location.start)
            end = int(genes.location.end)
            collector_dict.update({id_index: [start, end]})
            try:
                gene_names.update({id_index:[genes.sub_features[0].qualifiers['gene'],
                                             genes.sub_features[0].qualifiers['product']]})
            finally:
                continue

# Remove irrelevant entries
collector_dict.pop('NC_018406.1:1..964110', None)
collector_dict.pop('id-NC_018406.1:317779..317884', None)
collector_dict.pop('id-NC_018406.1:888405..890816', None)

# Initialize DataFrame to collect VCF data
collection_vcf = pd.DataFrame(columns=['Sample', 'Reference', 'Alternative', 'Position', 'Effect'])

# Process VCF files
os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/snpEff_results_virginia/')
for file in os.listdir():
    if file.endswith('.vcf'):
        vcf_reader = pysam.VariantFile(file, 'r')
        filename = str(file).replace('.bam.vcf.gz.vcf', '')
        list_refs = []
        list_alts = []
        list_pos = []
        list_effects = []
        for record in vcf_reader:
            try:
                eff_info = record.info['EFF']
                for element in eff_info:
                    list_refs.append(record.ref)
                    list_alts.append(record.alts[0])
                    list_pos.append(record.pos)
                    list_effects.append(element)
            except KeyError:
                continue

        adding_df = pd.DataFrame({
            'Sample': filename,
            'Reference': list_refs,
            'Alternative': list_alts,
            'Position': list_pos,
            'Effect': list_effects
        })
        collection_vcf = pd.concat([collection_vcf, adding_df], ignore_index=True)

# Define a function to find the corresponding gene ID for SNP position
def find_gene_id(position, gene_dict):
    for gene_id, (beginning, end_pos) in gene_dict.items():
        if beginning <= position <= end_pos:
            return gene_id
    return None

# Function to map values using the dictionary, keeping original if not found
def map_values(value, mapping):
    return mapping.get(value, value)

# Apply the function to add the 'Gene_ID' column
collection_vcf['Gene_ID'] = collection_vcf['Position'].apply(lambda pos: find_gene_id(pos, collector_dict))
collection_vcf['Gene_product'] = collection_vcf['Gene_ID'].apply(lambda x: map_values(x, gene_names))

# Filter non-synonymous SNPs
non_synonymous = collection_vcf[collection_vcf['Effect'].str.contains('STOP') ]# |
                               # (collection_vcf['Effect'].str.contains('NON_SYNONYMOUS'))]

# Define strains
motility_strains = ['A1', 'F1', 'F4', 'A10', 'E11', 'E12']
non_motility_strains = ['B2', 'A9', 'D8', 'C3', 'B8']

# Separate non-synonymous SNPs into motility and non-motility strains
non_synonymous_motility = non_synonymous[non_synonymous['Sample'].isin(motility_strains)]
non_synonymous_non_motility = non_synonymous[non_synonymous['Sample'].isin(non_motility_strains)]

# Identify genes present in non-motility strains
genes_in_non_motility = set(non_synonymous_non_motility['Gene_ID'])

# Identify genes present in motility strains
genes_in_motility = set(non_synonymous_motility['Gene_ID'])

# Find genes that are in non-motility strains but not in motility strains
genes_only_in_non_motility = genes_in_non_motility.difference(genes_in_motility)

# Filter the DataFrame to keep only these gene IDs
filtered_df = non_synonymous_non_motility[non_synonymous_non_motility['Gene_ID'].isin(genes_only_in_non_motility)]

# List of gene IDs that are non-synonymous and in non-motile strains but not in motile strains
list_genes_non_motility = list(sorted(set(filtered_df['Gene_ID'])))

# %% Creation of fasta file of genes of interest.
with open('Non_motile_variant_proteins_stop_codon.fasta', 'a') as handle:
    for genome in SeqIO.parse('/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCF_000286675.1/genomic.gbff', 'genbank'):
        for feature in genome.features:
            if feature.type == "CDS":  # Find CDS to collect the information
                locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '_')
                if locustag in list_genes_non_motility:
                    if "translation" in feature.qualifiers.keys():
                        aaseq = Seq(feature.qualifiers["translation"][0])
                        fasta_aa = SeqRecord(aaseq, locustag, description='')
                        SeqIO.write(fasta_aa, handle, 'fasta')

#%% Create student matrix
student_matrix = pd.DataFrame(columns=list(collection_vcf['Sample'].unique()), index=list(collector_dict.keys())).fillna(0)

for index, row in collection_vcf.iterrows():
    try:
        student_matrix[row['Sample']][row['Gene_ID']] += 1
    finally:
        continue

student_matrix.to_csv('Matrix_mutations_SNPs.tsv', sep='\t')