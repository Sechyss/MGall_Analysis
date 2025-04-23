#%% Load libraries and data
import pickle
import pandas as pd


base_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_HF'
presence_absence_file = pd.read_csv(f'{base_path}/gene_presence_absence_filt_pseudo_length_frag.csv', sep=',', header=0)

#%% Lipoprotein HMM scan
hmm_table = pd.read_csv(f'{base_path}/Lineage_differences/Lipoprotein/HMM_scan/lipoprotein_scan_seqs_simple.tsv', sep='\t', header=0)

lipoproteins = hmm_table['query_name'].unique()

filtered_presence_absence = presence_absence_file[
    presence_absence_file.apply(lambda row: any(item in lipoproteins for item in row), axis=1)
]

cluster_lipoproteins = filtered_presence_absence['Gene'].unique()
# Create a dictionary to store the clusters and their corresponding lipoproteins
cluster_dict = {'Lipoprotein' : cluster_lipoproteins}
with open(f'{base_path}/Lineage_differences/Lipoprotein/HMM_scan/lipoprotein_clusters.pickle', 'wb') as handle:
    pickle.dump(cluster_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
# %%

base_path = '/home/albertotr/OneDrive/Data/Ipoutcha_motility/BLAST_db'
blastp_db = pd.read_table(f'{base_path}/candidates_motility_pangenome.tsv', sep='\t', header=None)

motility = blastp_db[1].unique()
filtered_presence_absence = presence_absence_file[
    presence_absence_file.apply(lambda row: any(item in motility for item in row), axis=1)
]
cluster_motility = filtered_presence_absence['Gene'].unique()
# Create a dictionary to store the clusters and their corresponding motility genes
cluster_dict.update({'Motility' : cluster_motility})
with open(f'{base_path}/candidates_motility_pangenome.pickle', 'wb') as handle:
    pickle.dump(cluster_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
# %%
