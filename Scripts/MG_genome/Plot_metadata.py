from gettext import dpgettext
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/')

# Style settings (keep context but refine aesthetics)
sns.set_theme(style="whitegrid", context="talk")
sns.set_palette("Set2")
plt.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": 220,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "axes.titlesize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

metadata_GM = pd.read_excel('Metadata_genomes.xlsx', sheet_name='Metadata_Lucy', engine='openpyxl')
metadata_NCBI = pd.read_excel('Metadata_genomes.xlsx', sheet_name='Metadata_NCBI', engine='openpyxl')
metadata_NCBI = metadata_NCBI[(metadata_NCBI['Host'] == 'House finch') & (metadata_NCBI['index'].str.contains('GCA'))]
metadata_SRA = pd.read_excel('Metadata_genomes.xlsx', sheet_name='Metadata_SRA', engine='openpyxl')
metadata_SRA = metadata_SRA[metadata_SRA['Host'] == 'Haemorhous mexicanus']

# Keep your original arrays intact
dates = np.concatenate([metadata_GM['year.sampling'].values,
                        metadata_NCBI['Date'].values,
                        metadata_SRA['Collection_Date'].values])

location_origins = np.concatenate([metadata_GM['site.of.sampling'].values,
                                  metadata_NCBI['State'].values,
                                  metadata_SRA['geo_loc_name'].values])

def clean_locations(s: pd.Series) -> pd.Series:
    out = s.astype(str).str.strip()
    out = out.replace({'': np.nan, 'nan': np.nan, 'NA': np.nan, 'None': np.nan})
    out = out.dropna()
    out = out.str.title()
    return out

loc_parts = []
if 'site.of.sampling' in metadata_GM:
    loc_parts.append(clean_locations(metadata_GM['site.of.sampling']))
if 'State' in metadata_NCBI:
    loc_parts.append(clean_locations(metadata_NCBI['State']))
if 'geo_loc_name' in metadata_SRA:
    loc_parts.append(clean_locations(metadata_SRA['geo_loc_name']))

locations_all = pd.concat(loc_parts, ignore_index=True) if loc_parts else pd.Series(dtype=str)

loc_counts = locations_all.value_counts()
top_n = 12
if len(loc_counts) > top_n:
    top = loc_counts.iloc[:top_n].copy()
    other_sum = loc_counts.iloc[top_n:].sum()
    if other_sum > 0:
        top.loc['Other'] = other_sum
    loc_counts = top

# Plot histogram (years) and pie chart (locations)
fig, (ax_hist, ax_pie) = plt.subplots(1, 2, figsize=(14, 6))

# Beautified histogram of years (no change to data handling)
if dates is not None and len(dates) > 0:
    years_arr = np.asarray(list(dates))
    y_min, y_max = years_arr.min(), years_arr.max()
    bins = np.arange(y_min, y_max + 2)
    ax_hist.hist(years_arr, bins=bins, edgecolor='white', linewidth=0.8,
                 color=sns.color_palette("Set2")[0])
    ax_hist.set_xlabel('Sampling Year')
    ax_hist.set_ylabel('Count')
    ax_hist.set_title('Sampling Years')
    ax_hist.set_xticks(np.arange(y_min, y_max + 1))
    ax_hist.tick_params(axis='x', rotation=45)
    ax_hist.grid(axis='y', linestyle=':', alpha=0.5)
else:
    ax_hist.text(0.5, 0.5, 'No valid dates', ha='center', va='center')
    ax_hist.set_axis_off()

# Beautified pie chart of locations
if not loc_counts.empty:
    # Palette sized to number of slices
    palette = sns.color_palette("Set2", n_colors=len(loc_counts))
    # Slight explode for the largest slice to emphasize it
    largest_idx = int(np.argmax(loc_counts.values))
    explode = [0.05 if i == largest_idx else 0.02 for i in range(len(loc_counts))]

    # No values on the wedges; rely on legend only
    wedges = ax_pie.pie(
        loc_counts.values,
        labels=None,
        startangle=90,
        counterclock=False,
        explode=explode,
        colors=palette,
        wedgeprops=dict(edgecolor='white', linewidth=0.8)
    )[0]

    # Add a clean legend with counts
    legend = ax_pie.legend(
        wedges,
        [f"{lab} ({cnt})" for lab, cnt in zip(loc_counts.index, loc_counts.values)],
        title="Location",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=9,
        frameon=False
    )
    legend.get_title().set_fontsize(10)

    ax_pie.set_title('Location of Isolation', pad=12)
    ax_pie.axis('equal')
else:
    ax_pie.text(0.5, 0.5, 'No valid locations', ha='center', va='center')
    ax_pie.set_axis_off()

plt.tight_layout()
plt.savefig('Metadata_genomes_plot.png', bbox_inches='tight', dpi=600)
plt.show()
