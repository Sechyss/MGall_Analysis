import pandas as pd

base_path = '/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_HF'
presence_absence_file = f'{base_path}/gene_presence_absence_filt_pseudo_length_frag.csv'

hmm_table = f'{base_path}/Lineage_differences/Lipoporotein/HMM_scan/lipoprotein_scan_seqs_simple.tsv'
hmm_table = pd.read_csv(hmm_table, sep='\t', header=0)