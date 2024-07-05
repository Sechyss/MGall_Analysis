from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

os.chdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/pangenome_results_filtered/aligned_gene_sequences/')
for file in os.listdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/pangenome_results_filtered/aligned_gene_sequences/'):
    fastafile = SeqIO.parse(file, 'fasta')
    with open(str(file).replace('.aln.fas', '.fasta').replace('~', '-'), 'a') as handle:
        for record in fastafile:
            header = str(record.id).split(';')[0]+'_'+str(file).replace('.aln.fas', '').replace('~', '-')
            sequence = record.seq
            seq_record = SeqRecord(sequence, id=header, description='')
            SeqIO.write(seq_record, handle, 'fasta')
