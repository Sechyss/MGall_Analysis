import os
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
from Bio.SeqRecord import SeqRecord

#%% Creation of the database of proteome

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_SPADES_denovo')

with open('/home/albertotr/OneDrive/Data/Ipoutcha_motility/Lucy_protein_db.fna', 'w') as outfile:
    for directory in tqdm(os.listdir()):
        if directory.endswith('_prokka'):
            filename = str(directory)+'/'+ str(directory).replace('_prokka', '')+'.ffn'
            fasta_file = SeqIO.parse(filename, 'fasta')
            for record in fasta_file:
                sequence_id = record.id
                sequence = record.seq
                description = record.description
                fasta_aa = SeqRecord(sequence, sequence_id, description=description)
                SeqIO.write(fasta_aa, outfile, 'fasta')

#%% Filtering the database based on BLASTp
# Parse the FASTA file once and convert to list (this allows reuse)
fasta_db = list(SeqIO.parse('/home/albertotr/OneDrive/Data/Ipoutcha_motility/Lucy_protein_db.fna', 'fasta'))

blastp_table = pd.read_table('/home/albertotr/OneDrive/Data/Ipoutcha_motility/BLAST_db/BLASTP_Candidates_motility_newList.tsv', sep='\t', header=None)

# Initialize an empty dictionary
result_dict = {}

# Iterate through the DataFrame and populate the dictionary
for key, value in zip(blastp_table[0], blastp_table[1]):
    if key in result_dict:
        result_dict[key].append(value)
    else:
        result_dict[key] = [value]

#%% New Fasta file
for key in tqdm(result_dict.keys()):
    with open('/home/albertotr/OneDrive/Data/Ipoutcha_motility/BLAST_db/Homologous_'+ key +'.fna','a') as outfile:
        list_homologous = result_dict[key]
        for record in fasta_db:
            if record.id in list_homologous:
                fasta_aa = SeqRecord(record.seq, record.id, description=record.description)
                SeqIO.write(fasta_aa, outfile, 'fasta')
    outfile.close()
