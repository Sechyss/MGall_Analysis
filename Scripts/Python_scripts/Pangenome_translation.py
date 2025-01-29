from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

pangenome = SeqIO.parse("/home/albertotr/OneDrive/Data/Cambridge_Project/"
                        "pangenome_results_filtered/pan_genome_reference.fa", "fasta")

with open('/home/albertotr/Analysis/HMM_analysis/Lineage_differences/pangenome_reference.faa', 'a') as f1:
    for gene in pangenome:
        sequence = gene.seq
        prot_sequence = sequence.translate(table=4, stop_symbol='')
        seq_record_2 = SeqRecord(prot_sequence, id=gene.id, description='')
        SeqIO.write(seq_record_2, f1, 'fasta')
