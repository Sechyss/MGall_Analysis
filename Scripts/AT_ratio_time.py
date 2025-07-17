from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# === Parameters ===
alignment_file = "alignment.fasta"
metadata_file = "metadata.csv"
output_plot = "at_ratio_vs_time.png"

# === Load metadata ===
metadata = pd.read_csv(metadata_file)
metadata['strain'] = metadata['strain'].astype(str)

# === Compute AT ratios from alignment ===
at_data = []

for record in SeqIO.parse(alignment_file, "fasta"):
    strain_id = record.id
    seq = str(record.seq).upper().replace("-", "")  # remove gaps if any
    if len(seq) == 0:
        continue
    a = seq.count('A')
    t = seq.count('T')
    g = seq.count('G')
    c = seq.count('C')
    total = a + t + g + c
    at_ratio = (a + t) / total if total > 0 else None

    at_data.append({'strain': strain_id, 'at_ratio': at_ratio})

at_df = pd.DataFrame(at_data)

# === Merge AT ratios with metadata ===
df = pd.merge(metadata, at_df, on='strain')

# === Plot ===
plt.figure(figsize=(10, 6))
sns.scatterplot(data=df, x='date', y='at_ratio', s=60)
sns.regplot(data=df, x='date', y='at_ratio', scatter=False, color='red', ci=None)
plt.xlabel("Collection Year")
plt.ylabel("AT Ratio (A+T / Total)")
plt.title("AT Ratio Change Over Time")
plt.tight_layout()
plt.savefig(output_plot, dpi=300)
plt.show()
