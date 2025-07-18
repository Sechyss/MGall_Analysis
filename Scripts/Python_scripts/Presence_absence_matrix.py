#%% import of packages and function
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from ete3 import Tree


def filter_presence_absence(dataframe, group1, group2, filter1, filter2):
    dataframe_group1 = [col for col in dataframe.columns if col in group1]
    dataframe_group2 = [col for col in dataframe.columns if col in group2]

    # Calculate the presence ratio for house finch and poultry
    presence_group1 = dataframe[dataframe_group1].notna().mean(axis=1)
    presence_group2 = dataframe[dataframe_group2].notna().mean(axis=1)

    # Filter genes present in at least 90% of house finch samples and less than 10% of poultry samples
    filtered_genes_group1 = dataframe[
        (presence_group1 >= filter1) & (presence_group2 <= filter2)]
    filtered_genes_group2 = dataframe[
        (presence_group1 <= filter2) & (presence_group2 >= filter1)]
    return filtered_genes_group1, filtered_genes_group2

def fastafile_creation(fasta_file, list_genes, reference_file):
    pangenome_reference = SeqIO.parse(reference_file, 'fasta')
    with open(fasta_file+'.fasta', 'a') as f2:
        with open(fasta_file+'.faa', 'a') as f1:
            for gene in pangenome_reference:
                if gene.id in list_genes:
                    sequence = gene.seq
                    prot_sequence = sequence.translate(table=4, stop_symbol='')
                    seq_record = SeqRecord(sequence, id=gene.id, description='')
                    seq_record_2 = SeqRecord(prot_sequence, id=gene.id, description='')
                    SeqIO.write(seq_record, f2, 'fasta')
                    SeqIO.write(seq_record_2, f1, 'fasta')

Presence_absence = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                               'pangenome_results_HF/gene_presence_absence_filt_pseudo_length_frag.csv')

#%% Analysis presence/absence
lineage2 = [
'A090809_2009',
'G_2015',
'L_2015',
'J_2015',
'A01546_2007',
'A01505_2007',
'A01510_2007',
'A565_2002',
'A2261_2000',
'A506_2002',
'A454_2001',
'A442_2001',
'A471_2001',
'A267_2000',
'A3611_2001',
'A509_2002',
'A277_2001',
'BP421_2002',
'A933_2001',
'A297_2001',
'A371_2001',
'A877_2003',
'MG31_AL_11_2011',
'A020_2011',
'OY79_2013',
'WU47_2013',
'A_2015',
'Q_2015',
'Black_2012',
'GW77_2013',
'E054_2013',
'KB64_2013',
'A1_2014',
'MG27_AL_11_2011',
'A046_2011',
'A3278_2012',
'A3240_2012',
'A3244_2012',
'A3225_2012',
'A2486_2012',
'A304_2003',
'A195_2003',
'A012_2011',
'MG28_AL_11_2011',
'A044_2011',
'A001_2011',
'MG22_AL_11_2011',
'A004_2011',
'A3316_2012',
'Blue_2012',
'A021_2011',
'MG23_AL_11_2011',
'A006_2011',
'MG25_AL_11_2011',
'MG32_AL_11_2011',
'AMGO_2011',
'MG29_AL_11_2011',
'A013_2011',
'MG26_AL_11_2011',
'A018_2011',
'MG30_AL_11_2011',
'NC08_2008.031-4-3P',
'NC06_2006.080-5-2P'
            ]

tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold_50_combined.finaltree.newick')
leaves = tree_file.get_leaves()

lineage1 = [leaf.name.replace("'", "") for leaf in leaves if leaf.name.replace("'", "") not in lineage2]
lineage1 = [elem.replace('_2011', '') if 'MG' in elem and '_AL_11_2011' in elem else elem for elem in lineage1]
lineage2 = [elem.replace('_2011', '') if 'MG' in elem and '_AL_11_2011' in elem else elem for elem in lineage2]


#%% filtering by lineage

lineage1_genes, lineage2_genes = filter_presence_absence(Presence_absence, lineage1, lineage2, 0.15, 0.0)

list_of_genes = lineage1_genes['Gene'].to_list()

pangenome_file = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/pan_genome_reference.fa'
output_fasta = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/Lineage_1_shell_differences'
fastafile_creation(output_fasta, list_of_genes, pangenome_file)