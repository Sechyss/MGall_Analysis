"""
Compile and summarise viral/defence-system annotations across all samples.

Reads PADLOC (anti-viral defence) and VIBRANT (prophage) output files for
all sequenced house-finch Mycoplasma gallisepticum genomes and merges the
results with sample metadata.  Produces consolidated tables of defence
system and prophage presence/absence across isolates.

Usage:
    Update the base path and file paths inside the script, then run:
        python Metadata_Viral_analysis.py

Outputs:
    - Merged DataFrames of defence system and prophage annotations per sample
      (written to Excel or CSV as configured at the bottom of the script)
"""

import pandas as pd
import os
import glob
import pickle
import warnings

# Suppress FutureWarning messages
warnings.simplefilter(action='ignore', category=FutureWarning)

key_data = pd.read_excel('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Metadata_genomes.xlsx',
                         sheet_name='Metadata_Keys')
key_data['Sample Name'] = key_data['Sample Name'].replace(
    {'/': '_', ' ': '_', '\\.': '_', '-': '_', '\(': '_', '\)': ''}
    , regex=True)
key_data.set_index('Sample Name', inplace=True)

#%% Analysis of padloc for antiviral mechanisms

columns_table = ['system.number', 'seqid', 'system', 'target.name', 'hmm.accession',
                 'hmm.name', 'protein.name', 'full.seq.E.value', 'domain.iE.value',
                 'target.coverage', 'hmm.coverage', 'start', 'end', 'strand',
                 'target.description', 'relative.position', 'contig.end', 'all.domains',
                 'best.hits']
collector_data = pd.DataFrame(columns=columns_table)
os.chdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Padloc_results/')
for file in os.listdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Padloc_results/'):
    if file.endswith('.csv'):
        filename = file.replace('_padloc.csv', '')
        table = pd.read_csv(file)
        table['Genome'] = str(filename)
        table['Host'] = key_data.loc[filename, 'Host']
        table['Date'] = key_data.loc[filename, 'Date']
        collector_data = pd.concat([collector_data, table], ignore_index=True)

filtering_dict = pickle.load(
    open('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/filtering_dict.pickle', 'rb'))

#%% Extraction of data from VIBRANT

collector_data_2 = pd.DataFrame()
for file in glob.glob('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/'
                      'CheckMbins/*_vibrant/VIBRANT_*/VIBRANT_results_*/VIBRANT_summary_results_*tsv'):
    table = pd.read_table(file, sep='\t')
    filename = str(os.path.basename(file)).replace('VIBRANT_summary_results_', '').replace('.tsv', '')
    table['Genome'] = str(filename)
    table['Host'] = key_data.loc[filename, 'Host']
    table['Date'] = key_data.loc[filename, 'Date']
    table_2 = pd.read_table(str(file).replace('VIBRANT_summary_results_', 'VIBRANT_genome_quality_'), sep='\t')
    new_df = table.merge(table_2, how='outer', on='scaffold')
    collector_data_2 = pd.concat([collector_data_2, new_df], ignore_index=True)

writer = pd.ExcelWriter('/Users/at991/OneDrive - University of Exeter/Data/'
                        'Cambridge_Project/Metadata_phage_analysis.xlsx', engine='openpyxl')
collector_data.to_excel(writer, index=False, sheet_name='Padloc_analysis')
collector_data_2.to_excel(writer, index=False, sheet_name='Vibrant_analysis')
writer.close()
