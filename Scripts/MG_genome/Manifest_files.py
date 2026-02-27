"""
Generate NCBI submission manifest files with MD5 checksums.

Scans the raw-reads directory for 'Sample_*' folders containing paired-end
FASTQ files, translates folder names to standardised sample identifiers using
the Camille replacement dictionary, and builds a manifest DataFrame in the
format required for ENA/SRA submission.  A second step reads pre-computed
MD5 sums from an Excel table and incorporates them into the manifest.

Usage:
    Update 'base' and associated paths, then run:
        python Manifest_files.py

Outputs:
    - MGall_Manifest.csv           Initial manifest with filename placeholders
    - MGall_Manifest_with_MD5.csv  Manifest with actual MD5 checksums
"""

#%% Import packages
import os
import pandas as pd
import pickle

base = '/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/Raw_reads/'

#%%
with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements_foldername.pickle', 'rb') as handle:
    replacements = pickle.load(handle)

dictionary_manifest = {}
for folder in os.listdir(base):
    if folder.startswith('Sample_'):
        name = folder.replace('Sample_', '')
        new_name = replacements.get(name, name)
        folder_path = os.path.join(base, folder)
        fastq_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.fastq')])
        # Assume paired-end: first is forward, second is reverse
        if len(fastq_files) >= 2:
            forward_file, reverse_file = fastq_files[:2]
        else:
            forward_file = fastq_files[0] if fastq_files else ''
            reverse_file = ''
        dictionary_manifest[new_name] = [
            'Mgall_GNM_2025', 'Illumina HiSeq 1000', '', 'GENOMIC', 'unspecified', 'WGS', 'PAIRED',
            forward_file, forward_file+'.md5', reverse_file, reverse_file+'.md5'
        ]

manifest_df = pd.DataFrame.from_dict(dictionary_manifest, orient='index')
manifest_df.columns = [
    'study', 'instrument_model', 'library_name', 'library_source', 'library_selection',
    'library_strategy', 'library_layout', 'forward_file_name', 'forward_file_md5',
    'reverse_file_name', 'reverse_file_md5'
]

manifest_df.index.name = 'sample'
manifest_df.reset_index(inplace=True)
manifest_df.to_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/MGall_Manifest.csv', index=False)

#%% Replace MD5 with actual MD5 sums
old_df = pd.read_table('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/Raw_reads/fastq2_template_1756715460962.tsv', skiprows=1)

#Dictionary
md5_values = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/md5_values.xlsx', header=0, sheet_name='Sheet1')

# Create MD5 dictionary: {filename: md5sum}
md5_dict = dict(zip(md5_values['MD5_file'], md5_values['value']))

# Replace MD5 file names in old_df with actual MD5 sums
for col in ['forward_file_md5', 'reverse_file_md5']:
    old_df[col] = old_df[col].map(md5_dict).fillna(old_df[col])

old_df.to_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/MGall_Manifest_with_MD5.csv', index=False)

# %%
