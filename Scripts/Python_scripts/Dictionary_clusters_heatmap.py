#%% Load libraries and data
import pickle
import os
import pandas as pd


base_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_HF'
presence_absence_file = pd.read_csv(f'{base_path}/gene_presence_absence_filt_pseudo_length_frag.csv', sep=',', header=0)

# Lipoprotein HMM scan
hmm_table = pd.read_csv(f'{base_path}/Lineage_differences/Lipoprotein/HMM_scan/lipoprotein_scan_seqs_simple.tsv', sep='\t', header=0)

lipoproteins = hmm_table['query_name'].unique()

filtered_presence_absence = presence_absence_file[
    presence_absence_file.apply(lambda row: any(item in lipoproteins for item in row), axis=1)
]

cluster_lipoproteins = filtered_presence_absence['Gene'].unique()
# Create a dictionary to store the clusters and their corresponding lipoproteins
cluster_dict = {'Lipoprotein' : list(cluster_lipoproteins)}

# Other virulence factors

cluster_virulence = []

for file in os.listdir(f'{base_path}/Lineage_differences/Virulence_genes/HMM_scans/'):
    if file.endswith('.tsv'):
        # Read the HMM scan file
        hmm_table = pd.read_table(f'{base_path}/Lineage_differences/Virulence_genes/HMM_scans/{file}', sep='\t', header=0)
        # Extract the query names
        query_names = hmm_table['query_name'].unique()
        # Filter the presence-absence file for the query names
        filtered_presence_absence = presence_absence_file[
            presence_absence_file.apply(lambda row: any(item in query_names for item in row), axis=1)
        ]
        # Get the unique genes for this virulence factor
        unique_genes = list(filtered_presence_absence['Gene'].unique())
        # Add to the dictionary with the file-specific key
        cluster_dict[file.split('_')[0]] = unique_genes
        # Append the unique genes to the cluster_virulence list
        cluster_virulence.append(unique_genes)

# Flatten the cluster_virulence list of lists
cluster_virulence = [gene for sublist in cluster_virulence for gene in sublist]

# Cas9 genes

cas9_genes = presence_absence_file[presence_absence_file['Gene'].str.contains('cas9', na=False)]
cluster_cas9 = cas9_genes['Gene'].unique()
# Add Cas9 genes to the dictionary
cluster_dict.update({'Cas9': list(cluster_cas9)})

# Motility genes

base_path_2 = '/home/albertotr/OneDrive/Data/Ipoutcha_motility/BLAST_db'
blastp_db = pd.read_table(f'{base_path_2}/candidates_motility_pangenome.tsv', sep='\t', header=None)

motility = blastp_db[1].unique()
filtered_presence_absence = presence_absence_file[
    presence_absence_file.apply(lambda row: any(item in motility for item in row), axis=1)
]
cluster_motility = filtered_presence_absence['Gene'].unique()
# Create a dictionary to store the clusters and their corresponding motility genes
cluster_dict.update({'Motility' : list(cluster_motility)})

# %% Check for any redundant genes in the clusters and select only unique ones

unique_to_lipoproteins = set(cluster_dict['Lipoprotein']).difference(cluster_dict['Motility'])
unique_to_motility = set(cluster_dict['Motility']).difference(cluster_dict['Lipoprotein'])

cluster_dict['Lipoprotein'] = list(unique_to_lipoproteins)
cluster_dict['Motility'] = list(unique_to_motility)

# %% Combine all the clusters into a single list
# Combine all the clusters into a single list
all_clusters = set(cluster_dict['Lipoprotein'] + cluster_dict['Motility'] + cluster_virulence + list(cluster_cas9))

# Update the cluster dictionary with the key 'Virulence'
cluster_dict['Virulence'] = list(all_clusters)

# Save the dictionary to a pickle file

with open(f'{base_path}/candidates_clusters_pangenome.pickle', 'wb') as handle:
    pickle.dump(cluster_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
# %%
