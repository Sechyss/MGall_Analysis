from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def sequence_dictionary(fastafile):
    yielding_dict = {}
    for record in SeqIO.parse(fastafile, "fasta"):
        sequence = record.seq
        header = record.id
        yielding_dict[header] = sequence
    return yielding_dict


fasta_masked_alb = sequence_dictionary('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                               'VCF_trees/Gubbins_Rlow_onlyLucy/Rlow_gubbins.masked.aln')
fasta_masked_lcy = sequence_dictionary('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                                       'VCF_trees/Rlow_Lucy.masked.aln')

for keys in tqdm(fasta_masked_lcy.keys()):
    with open('/home/albertotr/OneDrive/Data/Cambridge_Project/'
                                       'VCF_trees/Lucy_alb_alignment_comparison/'+keys+'.fasta','a') as fasta_out:
        alb_id = str(keys)+'_alb'
        lcy_id = str(keys)+'_lcy'
        sequence_alb = fasta_masked_alb[keys]
        sequence_lcy = fasta_masked_lcy[keys]
        record_1 = SeqRecord(sequence_alb, id=alb_id, description='')
        record_2 = SeqRecord(sequence_lcy, id=lcy_id, description='')
        SeqIO.write(record_1, fasta_out, 'fasta')
        SeqIO.write(record_2, fasta_out, 'fasta')
