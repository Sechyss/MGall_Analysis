from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=True, type=str)
    parser.add_argument('--output', required=True, type=str)
    params = parser.parse_args()

    pangenome = SeqIO.parse(params.input, "fasta")

    with open(params.output, 'a') as f1:
        for gene in pangenome:
            sequence = gene.seq
            prot_sequence = sequence.translate(table=4, stop_symbol='')
            seq_record_2 = SeqRecord(prot_sequence, id=gene.id, description='')
            SeqIO.write(seq_record_2, f1, 'fasta')
