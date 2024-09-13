import os

import pandas as pd
import pysam
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

collector_dict = {}
gene_names = {}
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


collector_dict.pop('NC_018406.1:1..964110')
collection_vcf = pd.DataFrame(columns=['Sample', 'Reference', 'Alternative', 'Position', 'Effect'])
os.chdir('/home/albertotr/OneDrive/Data/'
                               'Cambridge_Project/Mapped_output_VA94_7994_1_7P/snpEff_results_virginia/')
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
            finally:
                continue

        adding_df = pd.DataFrame({'Sample': filename, 'Reference': list_refs, 'Alternative': list_alts,
                                                'Position': list_pos, 'Effect': list_effects})
        collection_vcf = pd.concat([collection_vcf, adding_df], ignore_index=True)

# Define a function to find the corresponding gene ID for a SNP position
def find_gene_id(position, gene_dict):
    for gene_id, (beginning, end_pos) in gene_dict.items():
        if beginning <= position <= end_pos:
            return gene_id
    return None

# Function to map values using the dictionary, keeping original if not found
def map_values(value, mapping):
    return mapping.get(value, value)



# Apply the function to add the 'Gene_ID' column to the non_synonymous DataFrame
collection_vcf['Gene_ID'] = collection_vcf['Position'].apply(lambda pos: find_gene_id(pos, collector_dict))
collection_vcf['Gene_product'] = collection_vcf['Gene_ID'].apply(lambda x: map_values(x, gene_names))

non_synonymous = collection_vcf[(collection_vcf['Effect'].str.contains('STOP')) |
                                (collection_vcf['Effect'].str.contains('NON_SYNONYMOUS'))]

# %% Analysis of the data collected based on the Ipoutcha paper

locustag_motility = ['HFMG94VAA_RS00385', 'HFMG94VAA_RS01105', 'HFMG94VAA_RS02755', 'HFMG94VAA_RS03655',
                     'HFMG94VAA_RS03725', 'HFMG94VAA_RS04385']
motility_strains, non_motility_strains = ['A1', 'F1', 'F4', 'A10', 'E11', 'E12'], ['B2', 'A9', 'D8', 'C3', 'B8']
non_synonymous_motility = non_synonymous[non_synonymous['Sample'].isin(motility_strains)]
non_synonymous_non_motility = non_synonymous[non_synonymous['Sample'].isin(non_motility_strains)]

# Group by gene_id and count the number of unique strains
gene_strain_count = non_synonymous_non_motility.groupby('Gene_ID')['Sample'].nunique()

# Filter gene_ids that are present in exactly 5 unique strains
genes_in_5_strains = gene_strain_count[gene_strain_count == len(non_motility_strains)].index

# Filter the original DataFrame to keep only these gene_ids
filtered_df = non_synonymous_non_motility[non_synonymous_non_motility['Gene_ID'].isin(genes_in_5_strains)]

list_genes_motility = list(set(filtered_df['Gene_ID']))

# %% Creation of fasta file of genes of interest.
with open('Non_motile_variant_proteins.fasta', 'a') as handle:
    for genome in SeqIO.parse('/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCF_000286675.1/genomic.gbff', 'genbank'):
        for feature in genome.features:
            if feature.type == "CDS":  # Find CDS to collect the information
                locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '_')
                if locustag in list_genes_motility:
                    if feature.qualifiers["translation"][0]:
                        aaseq = Seq(feature.qualifiers["translation"][0])
                        fasta_aa = SeqRecord(aaseq, locustag, description='')
                        SeqIO.write(fasta_aa, handle, 'fasta')
