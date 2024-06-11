from Bio.SeqIO import parse
from tqdm import tqdm
import os
import pandas as pd


os.chdir('/Users/at991\OneDrive - University of Exeter\Data\Mgall_NCBI/ncbi_dataset\data')
dictionary_collection = {}
for directory in tqdm(os.listdir()):
    if 'G' in directory:
        try:
            genbank = parse(directory+'/genomic.gbff', 'genbank')
            for record in genbank:
                accession_id = record.id
                date = record.annotations['date']
                BioProject = record.dbxrefs[0].replace('BioProject:', '')
                BioSample = record.dbxrefs[1].replace('BioSample:', '')
                #Assembly = record.dbxrefs[2].replace('Assembly:', '')
                if 'references' in record.annotations.keys():
                    Authors = record.annotations['references'][0].authors
                    Journal = record.annotations['references'][0].journal
                    Journal2 = record.annotations['references'][1].journal
                else:
                    Authors = 'No listed authors'
                    Journal = 'Unpublished'
                    Journal2 = ''

            dictionary_collection.update({accession_id:[date, BioProject, BioSample, Authors, Journal, Journal2]})
        finally:
            continue

dataframe = pd.DataFrame.from_dict(dictionary_collection, orient='index', columns=['date', 'BioProject', 'BioSample',
                                   'Authors', 'Journal', 'Journal2'])
dataframe.to_csv('metadata_Mgall_NCBI_genbanks.csv')
