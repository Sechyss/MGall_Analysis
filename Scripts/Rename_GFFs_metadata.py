import os
import pandas as pd

os.chdir('/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/GFFs_genomes')

metadata_sra = pd.read_excel('../SraRunTable_filtered.xlsx', sheet_name='Filtered_Accessions')

for index, row in metadata_sra.iterrows():
    oldname = row['Run']
    newname = row['Sample Name'].replace(' ', '_')
    os.rename(oldname + '.gff', newname + '.gff')

#%% Change the name of genome files

metadata_genomes = pd.read_table('/Users/at991/OneDrive - University of Exeter/Data/Mgall_NCBI/'
                                 'ncbi_dataset/data/data_summary.tsv')
counter = 1
for index, row in metadata_genomes.iterrows():
    oldname = row['Assembly Accession']
    newname = row['Organism Qualifier'].replace('strain: ', '').replace(' ', '_').replace('/', '_')
    if oldname+'.gff' in os.listdir():
        if newname+'.gff' in os.listdir():
            os.rename(oldname + '.gff', newname + '_counter_'+str(counter)+'.gff')
            counter += 1
        else:
            os.rename(oldname + '.gff', newname + '.gff')
            counter = 0
    else:
        continue
