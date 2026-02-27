"""
Generate input file lists for Gubbins from mapped-read FASTA consensus files.

Scans directories containing per-sample *.bam.fasta consensus files (Lucy
internal samples and SRA-derived samples) and writes tab-delimited list files
in the format expected by Gubbins:
    <sample_id> <TAB> <absolute_path_to_fasta>

SRA sample names are translated to standardised identifiers using project
metadata.

Usage:
    Update the directory paths inside the script, then run:
        python Creation_list_input_gubbins.py

Outputs:
    - List_files_HF_VA94_7994_1_7P.txt     Lucy samples mapped to VA94 reference
    - List_files_HF_Rlow.txt               Lucy samples mapped to Rlow reference
    - List_files_HF_SRA_VA94_7994_1_7P.txt SRA samples mapped to VA94 reference
    - List_files_HF_SRA_Rlow.txt           SRA samples mapped to Rlow reference
"""

import os
import pandas as pd


def list_input_creation(filedirectory, outputfile):
    os.chdir(filedirectory)

    with open(outputfile, 'a') as handle:
        for file in os.listdir():
            if file.endswith('.bam.fasta'):
                header = str(file).replace('.bam.fasta', '')
                handle.write(f"{header}\t{filedirectory+file}\n")
    handle.close()


def list_input_creation_sra(filedirectory, outputfile, exchange_dict):
    os.chdir(filedirectory)
    with open(outputfile, 'a') as handle:
        for file in os.listdir():
            if file.endswith('.bam.fasta'):
                if str(file).replace('.bam.fasta', '') in exchange_dict.keys():
                    key_dict = str(file).replace('.bam.fasta', '')
                    header = str(exchange_dict[key_dict]).replace(' ', '_')
                    handle.write(f"{header}\t{filedirectory+file}\n")
    handle.close()


#%% Rename files from Lucy mapped reads

list_input_creation('/home/albertotr/OneDrive/Data/'
                    'Cambridge_Project/Mapped_output_VA94_7994_1_7P/', 'List_files_HF_VA94_7994_1_7P.txt')
list_input_creation('/home/albertotr/OneDrive/Data/'
                    'Cambridge_Project/Mapped_output_Rlow/', 'List_files_HF_Rlow.txt')

#%% Rename of files from SRA

metadata = pd.read_excel('/home/albertotr/OneDrive/Data//Cambridge_Project/Metadata_genomes.xlsx',
                         sheet_name='Metadata_SRA', index_col=0)

metadata_2 = metadata[metadata['Host'] == 'Haemorhous mexicanus']

dict_names = metadata_2['Sample Name'].to_dict()

list_input_creation_sra('/home/albertotr/OneDrive/Data/'
                        'Cambridge_Project/Mapped_output_SRA_VA94_7994_1_7P/', 'List_files_HF_SRA_VA94_7994_1_7P.txt',
                        dict_names)
list_input_creation_sra('/home/albertotr/OneDrive/Data/'
                        'Cambridge_Project/Mapped_output_SRA_Rlow/', 'List_files_HF_SRA_Rlow.txt',
                        dict_names)
