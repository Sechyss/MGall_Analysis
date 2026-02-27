"""
Calculate AT content ratio over time for Mycoplasma gallisepticum genomes.

Computes the AT ratio (A+T)/(A+T+G+C) from FASTA genome sequences, merges with
sample collection-date metadata, and produces a scatter/regression plot saved to disk.

Usage:
    Update the path variables (base, fasta_dir, metadata_file, etc.) to point
    to your data, then run:
        python AT_ratio_overtime.py

Outputs:
    - at_ratio_vs_time.png  (scatter + OLS regression plot)
    - OLS regression summary printed to console
"""

import glob
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

base = '/home/albertotr/OneDrive/Data/Cambridge_Project/'
fasta_dir = f'{base}/GWAS/HF_fna/'  # Directory containing FASTA files
metadata_file = f'{base}/Edited_Camille_all_samples_dates.csv'
dictionary_path = f'{base}/Camille_replacements_foldername.pickle' # Dictionary for strain replacements
# Load the dictionary for strain replacements
with open(dictionary_path, 'rb') as handle:
    replacements = pickle.load(handle)
output_plot = f'{base}/at_ratio_vs_time.png'

# === Load metadata ===
metadata = pd.read_csv(metadata_file)
metadata['strain'] = metadata['strain'].astype(str)

# === Compute AT ratios from multiple FASTA files ===
at_data = []
for fasta_path in glob.glob(f'{fasta_dir}/*.fna'):
    strain_id = fasta_path.split('/')[-1].split('.')[0]  # Adjust if needed
    a = t = g = c = 0
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq).upper().replace("-", "")
        a += seq.count('A')
        t += seq.count('T')
        g += seq.count('G')
        c += seq.count('C')
    total = a + t + g + c
    at_ratio = (a + t) / total if total > 0 else None
    at_data.append({'strain': strain_id, 'at_ratio': at_ratio})

at_df = pd.DataFrame(at_data)
# Replace strain names in at_df using the replacements dictionary
at_df['strain'] = at_df['strain'].map(replacements).fillna(at_df['strain'])

# === Merge AT ratios with metadata ===
df = pd.merge(metadata, at_df, on='strain')

# === Plot ===
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='date', y='at_ratio', s=60)
sns.regplot(data=df, x='date', y='at_ratio', scatter=False, color='red', ci=None)
plt.xlabel("Sample Year")
plt.ylabel("AT Ratio")
plt.tight_layout()
plt.savefig(output_plot, dpi=600)
plt.show()

import statsmodels.api as sm

X = df['date']
y = df['at_ratio']
X = sm.add_constant(X)
model = sm.OLS(y, X).fit()
print(model.summary())