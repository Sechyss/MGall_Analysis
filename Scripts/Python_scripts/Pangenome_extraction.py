import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
# Load the CSV file into a DataFrame
# Make sure to specify the correct path to your CSV file
# and the correct column names for your data

database = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Pangenome_results_HF/gene_data.csv', header=0)

# Use tqdm with total parameter to show progress
with tqdm(total=len(database), desc="Processing rows") as pbar:
    for index, row in database.iterrows():
        # Extract the sequence from the row
        sequence = Seq(row['prot_sequence'])
        
        # Create a SeqRecord object
        seq_record = SeqRecord(sequence, id=row['annotation_id'], description='')
        
        # Save the SeqRecord object to a file
        filename = '/home/albertotr/OneDrive/Data/Cambridge_Project/Pangenome_results_HF/ALL_proteins_pangenome.fasta'
        with open(filename, "a") as output_handle:
            SeqIO.write(seq_record, output_handle, "fasta")
        
        # Update the progress bar
        pbar.update(1)