"""
Translate a Panaroo pangenome nucleotide FASTA to protein sequences.

Reads a nucleotide FASTA produced by Panaroo (pan_genome_reference.fa or
similar) and translates each record using the Mycoplasma/Spiroplasma genetic
code (codon table 4, where UGA encodes tryptophan rather than a stop).
The translated protein sequences are appended to the output FASTA file.

Usage:
    python Pangenome_translation.py --input <pangenome.fasta> --output <proteins.fasta>

Arguments:
    --input   Path to the input pangenome nucleotide FASTA.
    --output  Path for the output protein FASTA (appended if file exists).

Outputs:
    - <output>  Protein FASTA with one record per input nucleotide sequence
"""

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
