#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pyreadr
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

data = pd.read_table('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/Skygrid_data.csv', sep='\t', skiprows=1)

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
data_lineages = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/Lineage_through_time.csv')
# Invert the data in data_lineages
data_lineages_inverted = data_lineages.iloc[::-1].reset_index(drop=True)

# Combine the inverted data with the original data
combined_df = pd.concat([data_lineages_inverted, data], axis=1)

# Select the lineages to plot
lineages_to_plot = [col for col in combined_df.columns if 'Group' in col]

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
plt.xlim(left=1991, right=time.max())
plt.ylim(bottom=0)
plt.tight_layout()

plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/Skygrid_Lineages_throughtime_per.png', dpi=600)

# Show the plot
plt.show()

# %% Plot Re_gridded_hpd1 and Re_gridded_hpd2
plt.figure(figsize=(10, 6))

plt.plot(Re_gridded_hpd1_df['times1'], Re_gridded_hpd1_df['med'], color='C0', label='Lineage1')
plt.fill_between(Re_gridded_hpd1_df['times1'], Re_gridded_hpd1_df['lower'], Re_gridded_hpd1_df['upper'], color='C0', edgecolor='black', alpha=0.2, label='95% CI')
plt.plot(Re_gridded_hpd2_df['times2'], Re_gridded_hpd2_df['med'], color='C1', label='Lineage2')
plt.fill_between(Re_gridded_hpd2_df['times2'], Re_gridded_hpd2_df['lower'], Re_gridded_hpd2_df['upper'], color='C1', edgecolor='black', alpha=0.2, label='95% CI')

# Set the x-axis limits to cover the range of both datasets
combined_time = pd.concat([Re_gridded_hpd1_df['times1'], Re_gridded_hpd2_df['times2']])
plt.xlim(left=combined_time.min(), right=combined_time.max())
#plt.xlim(left=time.min(), right=time.max())

plt.ylim(bottom=0, top=30)
# Remove the top and right spines (the square around the plot)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Add labels and title
plt.xlabel('Time')
plt.ylabel('Re')
plt.legend()
plt.tight_layout()

plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Re_gridded_hpd.png', dpi=600)

# Show the plot
plt.show()


#%% Combine the previous figures into a single PNG with A being the first figure and B the second one

# Create a new figure with two subplots (A and B)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Plot the first figure (A) on the first subplot (ax1)
ax1.plot(time, median_log, color='black')
ax1.fill_between(time, lower_log, upper_log, color='none', edgecolor='black', alpha=0.2, hatch='//')

bottom = np.zeros_like(time)
for lineage in lineages_to_plot:
    ax1.fill_between(time, bottom, bottom + combined_df[lineage] * median_log / 100, alpha=0.5, label=lineage.split('_')[0].replace('1', ' A').replace('2', ' B'))
    bottom += combined_df[lineage] * median_log / 100

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_ylabel('Log$_{10}$(N$_{e}$τ)')
ax1.legend()
ax1.set_xlim(left=1991, right=combined_time.max())
ax1.set_ylim(bottom=0)
ax1.set_title('A', loc='left')

# Plot the second figure (B) on the second subplot (ax2)
ax2.plot(Re_gridded_hpd1_df['times1'], Re_gridded_hpd1_df['med'], color='C0', label='Lineage1')
ax2.fill_between(Re_gridded_hpd1_df['times1'], Re_gridded_hpd1_df['lower'], Re_gridded_hpd1_df['upper'], color='C0', edgecolor='C0', alpha=0.2)
ax2.plot(Re_gridded_hpd2_df['times2'], Re_gridded_hpd2_df['med'], color='C1', label='Lineage2')
ax2.fill_between(Re_gridded_hpd2_df['times2'], Re_gridded_hpd2_df['lower'], Re_gridded_hpd2_df['upper'], color='C1', edgecolor='C1', alpha=0.2)

ax2.set_xlim(left=1991, right=combined_time.max())
ax2.set_ylim(bottom=0, top=30)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.set_xlabel('Time')
ax2.set_ylabel('R$_{e}$')
ax2.set_title('B', loc='left')

# Link the x-axes of both subplots
ax1.get_shared_x_axes().join(ax1, ax2)

plt.tight_layout()
plt.savefig('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Combined_Figures.png', dpi=600)

plt.show()

# %%
