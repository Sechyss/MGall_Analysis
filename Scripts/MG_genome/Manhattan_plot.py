#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import os

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/GWAS/')

# Set the font size and style for better readability
matplotlib.rcParams.update({'font.size': 12})
plt.style.use('seaborn-v0_8-whitegrid')

# Read the data
df = pd.read_csv('swelling_manhattan.txt', sep='\t')

# Define colors for alternating chromosomes
colors = ['#1f77b4', '#ff7f0e']  # blue and orange

# Calculate Bonferroni correction threshold
# (typically 0.05 / number of tests)
alpha = 0.05
n_tests = len(df)
bonferroni_threshold = -np.log10(alpha / n_tests)
print(f"Bonferroni threshold: {bonferroni_threshold}")

# Create the plot
plt.figure(figsize=(14, 8))

# Since all points are on one chromosome, we'll just use one color
plt.scatter(df['BP'], df['LOG10P'], alpha=0.8, s=15, color=colors[0])

# Add Bonferroni line
plt.axhline(y=bonferroni_threshold, color='red', linestyle='--', 
           label=f'Bonferroni (p={alpha/n_tests:.2e})')

# Highlight significant points (above Bonferroni threshold) if any
significant = df[df['LOG10P'] > bonferroni_threshold]
if not significant.empty:
    plt.scatter(significant['BP'], significant['LOG10P'], color='red', s=30, 
               label=f'Significant SNPs (n={len(significant)})')

# Add labels and title
plt.xlabel('Position (bp)')
plt.ylabel('-log10(p-value)')
plt.title('Manhattan Plot for Swelling GWAS')

# Set y-axis to start at 0
y_max = max(df['LOG10P'].max() * 1.05, bonferroni_threshold * 1.1)
plt.ylim([0, y_max])

# Add legend
plt.legend(loc='upper left')

# Tighten layout
plt.tight_layout()

# Save figure
plt.savefig('swelling_manhattan_plot.png', dpi=600)
print("Plot saved as swelling_manhattan_plot.png")
plt.savefig('swelling_manhattan_plot.svg', dpi=600)
print("Plot saved as swelling_manhattan_plot.svg")
plt.savefig('swelling_manhattan_plot.pdf', dpi=600)
print("Plot saved as swelling_manhattan_plot.pdf")

# Display the plot
