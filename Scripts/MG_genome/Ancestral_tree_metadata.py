"""
Convert a lineage-membership dictionary to a PastML/BEAST metadata CSV.

Loads a pickle dictionary mapping lineage labels to lists of sample IDs
(formatted as <sample>_<year>), extracts the year from each identifier, and
writes a CSV with columns [Taxa, Lineage, Year] indexed by taxon name.
The output is used as metadata input for ancestral-state reconstruction
(PastML) or dated Bayesian phylogenetic analysis (BEAST).

Usage:
    Update the pickle path and output CSV path, then run:
        python Ancestral_tree_metadata.py

Outputs:
    - lineage_metadata.csv  Metadata table with columns: Lineage, Year (Taxa as index)
"""

import pickle
import pandas as pd


pickle_in = open('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/lineage1.pickle', 'br')
lineage_dict = pickle.load(pickle_in)
pickle_in.close()
data = []
for key, value in lineage_dict.items():
    for item in value:
        parts = item.rsplit('_', 1)
        if len(parts) == 2:
            data.append({'Lineage': key, 'Year': parts[1], 'Taxa': item})

df = pd.DataFrame(data)
df.set_index('Taxa', inplace=True)
df.to_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage_metadata.csv')