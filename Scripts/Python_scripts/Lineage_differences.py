import pickle
import pandas as pd
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from venny4py.venny4py import *


def fastafile_creation(fasta_file, list_genes, reference_file):
    with open(fasta_file, 'a') as handle:
        for genome in SeqIO.parse(reference_file, 'genbank'):
            for feature in genome.features:
                if feature.type == "CDS":  # Find CDS to collect the information
                    locustag = str(feature.qualifiers["locus_tag"][0]).replace(' ', '_')
                    if locustag in list_genes:
                        if "translation" in feature.qualifiers.keys():
                            aaseq = Seq(feature.qualifiers["translation"][0])
                            fasta_aa = SeqRecord(aaseq, locustag, description='')
                            SeqIO.write(fasta_aa, handle, 'fasta')

dictionary_file_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/mgall_60threshold_VA94/lineage_dict.pickle'
lineage_dict = pickle.load(open(dictionary_file_path, 'rb'))

mutation_df = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/All_mutations_matrix.csv', index_col=0)
lineage_1 = lineage_dict[1]
lineage_2 = lineage_dict[2]

lineage_1_snps = mutation_df[mutation_df['Sample'].isin(lineage_1)]
lineage_2_snps = mutation_df[mutation_df['Sample'].isin(lineage_2)]

unique_to_lineage1 = lineage_1_snps[~lineage_1_snps['Position'].isin(lineage_2_snps['Position'])]
unique_to_lineage2 = lineage_2_snps[~lineage_2_snps['Position'].isin(lineage_1_snps['Position'])]

genes_affected_l1 = set(lineage_1_snps['Gene_ID'].to_list())
genes_affected_l2 = set(lineage_2_snps['Gene_ID'].to_list())

unique_genes_lineage2 = genes_affected_l2.difference(genes_affected_l1)
filtered_df = lineage_2_snps[lineage_2_snps['Gene_ID'].isin(unique_genes_lineage2)]

#%% Figures

sets = {
    'Lineage1_genes': {x for x in genes_affected_l1 if not (isinstance(x, float) and math.isnan(x))},
    'Lineage2_genes': {x for x in genes_affected_l2 if not (isinstance(x, float) and math.isnan(x))},
}

venny4py(sets=sets, out='/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/Different_lineages_VA94_Genes_snps_Venn.png',dpi=600)
fastafile_creation('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/Different_lineages_VA94_Genes_snps.faa',
unique_genes_lineage2,
'/home/albertotr/OneDrive/Data/MGall_NCBI/ncbi_dataset/data/GCF_000286675.1/genomic.gbff')

#%% Lineages Rlow

dictionary_file_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/PopPUNK/mgall_60threshold/lineage_dict.pickle'
lineage_Rlow_dict = pickle.load(open(dictionary_file_path, 'rb'))
