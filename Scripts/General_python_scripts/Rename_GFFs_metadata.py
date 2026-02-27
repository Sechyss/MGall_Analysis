"""
Rename genome FASTA files to use standardised sample identifiers.

Reads metadata from an SRA run table and an NCBI data-summary TSV to map
accession-based filenames to human-readable strain names, then renames the
.fna files in place.  A counter suffix is appended when name collisions occur.

Usage:
    Update the directory path (os.chdir) and metadata file paths inside the
    script, then run:
        python Rename_GFFs_metadata.py

Notes:
    - Expects .fna files named by assembly accession (e.g. GCA_XXXXXXXX.fna)
      in the working directory.
    - Spaces and forward-slashes in strain names are converted to underscores.
"""

import os
import pandas as pd

os.chdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/CheckMbins/')

metadata_sra = pd.read_excel('../SraRunTable_filtered.xlsx', sheet_name='Filtered_Accessions')

for index, row in metadata_sra.iterrows():
    oldname = row['Run']
    newname = row['Sample Name'].replace(' ', '_')
    os.rename(oldname + '.fna', newname + '.fna')

#%% Change the name of genome files

metadata_genomes = pd.read_table('/Users/at991/OneDrive - University of Exeter/Data/Mgall_NCBI/'
                                 'ncbi_dataset/data/data_summary.tsv')
counter = 1
for index, row in metadata_genomes.iterrows():
    oldname = row['Assembly Accession']
    newname = row['Organism Qualifier'].replace('strain: ', '').replace(' ', '_').replace('/', '_')
    if oldname+'.fna' in os.listdir():
        if newname+'.fna' in os.listdir():
            os.rename(oldname + '.fna', newname + '_counter_'+str(counter)+'.fna')
            counter += 1
        else:
            os.rename(oldname + '.fna', newname + '.fna')
            counter = 0
    else:
        continue
