import os
import pandas as pd
import pickle

base = '/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/Raw_reads/'

with open('/home/albertotr/OneDrive/Data/Cambridge_Project/Camille_replacements_foldername.pickle', 'rb') as handle:
    replacements = pickle.load(handle)

dictionary_manifest = {}
for folder in os.listdir(base):
    if folder.startswith('Sample_'):
        name = folder.replace('Sample_', '')
        new_name = replacements.get(name, name)  # Use .get() to avoid KeyError
        folder_path = os.path.join(base, folder)
        fastq_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.fastq.gz')])
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