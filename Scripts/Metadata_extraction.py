from Bio.SeqIO import parse
from tqdm import tqdm
import os
import pandas as pd


os.chdir('/Users/at991\OneDrive - University of Exeter\Data\Mgall_NCBI/ncbi_dataset\data')
data_summary = pd.read_table('data_summary.tsv', sep='\t', index_col='Assembly Accession')
dictionary_collection = {}
for directory in tqdm(os.listdir()):
    if 'GC' in directory:
        try:
            genbank = parse(directory+'/genomic.gbff', 'genbank')
            for record in genbank:
                accession_id = record.id
                date = record.annotations['date']
                if 'references' in record.annotations.keys():
                    Authors = record.annotations['references'][0].authors
                    Journal = record.annotations['references'][0].journal
                    Journal2 = record.annotations['references'][1].journal
                else:
                    Authors = 'No listed authors'
                    Journal = 'Unpublished'
                    Journal2 = ''

            dictionary_collection.update({directory:[date, accession_id, Authors, Journal, Journal2]})
        finally:
            continue

dataframe = pd.DataFrame.from_dict(dictionary_collection, orient='index', columns=['date', 'Locus_id',  'Authors',
                                                                                   'Journal',
                                                                                   'Journal2'])

new_df = dataframe.join(data_summary, how='outer')

new_df.to_csv('metadata_Mgall_NCBI_genbanks.csv')
