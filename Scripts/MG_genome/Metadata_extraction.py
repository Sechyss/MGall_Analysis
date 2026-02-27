"""
Extract and consolidate metadata from NCBI GenBank files and SRA run tables.

Parses GenBank flat files to retrieve accession, collection date, authors,
journal references, host, location, strain, and biosample IDs for each
sequenced Mycoplasma gallisepticum genome.  The extracted metadata is then
combined with the NCBI data-summary TSV, SRA run table, and internal project
keys and written to a multi-sheet Excel workbook.

Usage:
    Update os.chdir() and file paths inside the script, then run:
        python Metadata_extraction.py

Outputs:
    - Metadata_genomes.xlsx  Multi-sheet Excel workbook containing:
        * Metadata_NCBI   — GenBank-derived genome metadata
        * Metadata_SRA    — SRA run information
        * Metadata_Keys   — Merged metadata with internal project identifiers
"""

from Bio.SeqIO import parse
from tqdm import tqdm
import os
import pandas as pd

#%% Metadata from genbank files (NCBI)
os.chdir('/Users/at991/OneDrive - University of Exeter/Data/Mgall_NCBI/ncbi_dataset/data')
data_summary = pd.read_table('data_summary.tsv', sep='\t', index_col='Assembly Accession')
dictionary_collection = {}
for directory in tqdm(os.listdir()):
    if 'GC' in directory:
        try:
            genbank = parse(directory + '/genomic.gbff', 'genbank')
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

            dictionary_collection.update({directory: [date, accession_id, Authors, Journal, Journal2]})
        finally:
            continue

dataframe = pd.DataFrame.from_dict(dictionary_collection, orient='index', columns=['date', 'Locus_id', 'Authors',
                                                                                   'Journal',
                                                                                   'Journal2'])

new_df = dataframe.join(data_summary, how='outer').reset_index()

# new_df.to_csv('metadata_Mgall_NCBI_genbanks.csv')

#%% Combination of metadata with reads metadata

os.chdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project')
sra_table = pd.read_excel(
    '/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/SraRunTable_filtered.xlsx',
    sheet_name='Filtered_Accessions')

list_columns = ['Sample Name', 'Run', 'Collection_Date',  'Host']

sra_filtered = sra_table[list_columns]
sra_filtered.set_index('Sample Name', inplace=True)
sra_filtered.rename(columns={'Collection_Date': 'Date', 'Run': 'Locus_id'}, inplace=True)


new_df.replace({'strain: ': ''}, regex=True, inplace=True)
new_df2 = new_df[['Organism Qualifier', 'Locus_id', 'Submission Date', 'Host']].copy().reset_index()
new_df2.rename(columns={'Submission Date': 'Date', 'Organism Qualifier': 'Sample Name'}, inplace=True)
new_df2 = new_df2.sort_values(by=['Sample Name']).drop_duplicates(subset='Sample Name', keep='first')
new_df2['Date'] = pd.to_datetime(new_df2['Date'], format='%d/%m/%Y', dayfirst=True)
new_df2['Date'] = new_df2['Date'].dt.year
new_df2.drop(['index'], axis=1, inplace=True)
new_df2 = new_df2[new_df2['Sample Name'] != 'KUVMG001']
new_df2.set_index('Sample Name', inplace=True)

#%% Add Lucy's metadata

meta_table = pd.read_csv(
    '/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/meta.csv')

meta_filtered = meta_table[['trelabels', 'Folder.name', 'year.sampling']].copy()
meta_filtered.replace({'Sample_': ''}, regex=True, inplace=True)
meta_filtered.rename(columns={'trelabels': 'Sample Name', 'Folder.name': 'Locus_id', 'year.sampling': 'Date'},
                     inplace=True)
meta_filtered.set_index('Sample Name', inplace=True)
meta_filtered['Host'] = 'Haemorhous mexicanus'

final_df = pd.concat([sra_filtered, new_df2, meta_filtered], axis=0)

#%% Combination of all datasets into one single Excel spreadsheet

writer = pd.ExcelWriter('Metadata_genomes.xlsx', engine='openpyxl')
new_df.to_excel(writer, index=False, sheet_name='Metadata_NCBI')
sra_table.to_excel(writer, index=False, sheet_name='Metadata_SRA')
meta_table.to_excel(writer, index=False, sheet_name='Metadata_Lucy')
final_df.to_excel(writer, index=True, sheet_name='Metadata_Keys')
writer.close()
