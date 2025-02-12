#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pickle
from collections import defaultdict

data = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Skygrid_reconstruction.csv')

# Extract relevant columns
time = data['time']
median = data['median']
lower = data['lower']
upper = data['upper']
# Log the values
median_log = np.log10(median)
lower_log = np.log10(lower)
upper_log = np.log10(upper)

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(time, median_log, label='Median', color='blue')
plt.fill_between(time, lower_log, upper_log, color='grey', alpha=0.2, label='95% CI')

# Add labels and title
plt.xlabel('Time')
plt.ylabel('Log$_{10}$(N$_{e}$τ)')
plt.title('SkyGrid Plot')
plt.legend()

# Adjust the plot to remove space between the beginning of the lines and the y-axis
plt.xlim(left=time.min(), right=time.max())

# Show the plot
plt.show()

# %%
data_lineages = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Lineages_throughtime.csv')

pickle_in = open('/home/albertotr/OneDrive/Data/Cambridge_Project/pangenome_results_filtered/Lineage_differences/lineage1.pickle', 'rb')
lineages = pickle.load(pickle_in)

# Initialize a dictionary to hold the counts
lineage_counts = defaultdict(lambda: defaultdict(int))

# Iterate through the lineages dictionary
for lineage, taxa in lineages.items():
    for taxon in taxa:
        # Extract the date from the taxon string
        date_str = taxon.split('_')[-1]
        year = date_str.split('-')[0]
        # Increment the count for this lineage and year
        lineage_counts[lineage][year] += 1

# Create a set of all years
all_years = set()
for years in lineage_counts.values():
    all_years.update(years.keys())

# Create a dataframe
df = pd.DataFrame(index=sorted(all_years))

# Populate the dataframe with counts
for lineage, years in lineage_counts.items():
    df[lineage] = df.index.map(lambda year: years.get(year, 0))

# Reset the index to have 'Year' as a column
df.reset_index(inplace=True)
df.rename(columns={'index': 'Year'}, inplace=True)

print(df)
# %%
