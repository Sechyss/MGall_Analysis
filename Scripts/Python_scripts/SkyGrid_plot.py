#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyreadr

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

# Lineage differences over time and relative to the population
data_lineages = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Lineages_throughtime.csv')
# Invert the data in data_lineages
data_lineages_inverted = data_lineages.iloc[::-1].reset_index(drop=True)

# Combine the inverted data with the original data
combined_df = pd.concat([data_lineages_inverted, data], axis=1)

# Select the lineages to plot
lineages_to_plot = [col for col in combined_df.columns if 'percentage' in col]

# Open data from R in the format RData
results = pyreadr.read_r('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Re_data/lineages_data.RData')
Re_gridded_hpd1 = results['Re_gridded_hpd1'].transpose()
Re_gridded_hpd2 = results['Re_gridded_hpd2'].transpose()
times1 = results['times1']
times2 = results['times2']

Re_gridded_hpd1_df = pd.concat([times1.reset_index(drop=True), Re_gridded_hpd1.reset_index(drop=True)], axis=1)
Re_gridded_hpd2_df = pd.concat([times2.reset_index(drop=True), Re_gridded_hpd2.reset_index(drop=True)], axis=1)


# %% Estimation of the relative frequency of each lineage over time in the total population from data
# Create the plot
plt.figure(figsize=(10, 6))

# Plot the number of effective population over time
plt.plot(time, median_log, color='black')
plt.fill_between(time, lower_log, upper_log, color='none', edgecolor='black', alpha=0.2, hatch='//', label='95% CI')

# Fill the area under the curve using the relative frequency of the selected lineages
bottom = np.zeros_like(time)
for lineage in lineages_to_plot:
    plt.fill_between(time, bottom, bottom + combined_df[lineage] * median_log / 100, alpha=0.5, label=lineage)
    bottom += combined_df[lineage] * median_log / 100


# Remove the top and right spines (the square around the plot)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
# Add labels and title
plt.xlabel('Time')
plt.ylabel('Log$_{10}$(N$_{e}$τ)')
plt.legend()
# Adjust the plot to remove space between the beginning of the lines and the y-axis
plt.xlim(left=time.min(), right=time.max())
plt.ylim(bottom=0)
plt.tight_layout()

plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Skygrid_Lineages_throughtime_per.png', dpi=600)

# Show the plot
plt.show()

# %%

