#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
from matplotlib.ticker import ScalarFormatter
import os


def parse_gff(gff_file, target_chromosome=None):
    """
    Parse a GFF file and extract CDS features with gene information.
    
    Args:
        gff_file (str): Path to the GFF file
        target_chromosome (str, optional): If provided, only features from this chromosome are extracted
        
    Returns:
        list: List of dictionaries containing CDS features with their positions and annotations
    """
    cds_features = []
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
                
            chromosome, source, feature_type, start, end, score, strand, phase, attributes = parts
            
            # If target_chromosome specified, only keep features from that chromosome
            if target_chromosome and chromosome != target_chromosome:
                continue
                
            # We're only interested in CDS features
            if feature_type != 'CDS':
                continue
                
            # Extract gene name and product from attributes
            gene_name = "Unknown"
            product = "Unknown"
            
            for attribute in attributes.split(';'):
                if '=' in attribute:
                    key, value = attribute.split('=', 1)
                    if key == 'gene':
                        gene_name = value
                    elif key == 'product':
                        product = value
            
            cds_features.append({
                'chromosome': chromosome,
                'start': int(start),
                'end': int(end),
                'strand': strand,
                'gene_name': gene_name,
                'product': product
            })
    
    return cds_features


def find_nearby_genes(position, cds_features, distance=1000):
    """
    Find genes within specified distance of a position.
    
    Args:
        position (int): Position on the chromosome
        cds_features (list): List of CDS features
        distance (int): Maximum distance to consider a gene as nearby
        
    Returns:
        list: List of nearby gene features
    """
    nearby_genes = []
    
    for feature in cds_features:
        # Check if position is within feature or within distance of feature boundaries
        if (position >= feature['start'] - distance and 
            position <= feature['end'] + distance):
            
            # Calculate distance from position to feature
            if position >= feature['start'] and position <= feature['end']:
                dist = 0  # Position is within the feature
            else:
                dist = min(abs(position - feature['start']), abs(position - feature['end']))
                
            nearby_genes.append({
                'gene_name': feature['gene_name'],
                'product': feature['product'],
                'distance': dist,
                'strand': feature['strand']
            })
    
    # Sort by distance
    nearby_genes.sort(key=lambda x: x['distance'])
    
    return nearby_genes


def create_manhattan_plot(manhattan_data, gff_file, output_file='manhattan_plot.png', 
                          target_chromosome=None, distance=1000, pvalue_threshold=0.05):
    """
    Create a Manhattan plot with gene annotations for significant variants.
    
    Args:
        manhattan_data (str): Path to the Manhattan plot data file
        gff_file (str): Path to the GFF file
        output_file (str): Path to save the output plot
        target_chromosome (str, optional): If provided, only plot variants from this chromosome
        distance (int): Maximum distance to consider a gene as nearby
        pvalue_threshold (float): P-value threshold for significance (before Bonferroni correction)
    """
    # Read Manhattan plot data
    try:
        data = pd.read_csv(manhattan_data, sep='\t')
        print(f"Loaded {len(data)} variants from {manhattan_data}")
    except Exception as e:
        print(f"Error reading Manhattan data file: {e}")
        return
    
    # Ensure required columns are present
    required_columns = ['CHR', 'BP', 'P', 'LOG10P']
    if not all(col in data.columns for col in required_columns):
        print(f"Error: Manhattan data file must contain columns: {', '.join(required_columns)}")
        print(f"Found columns: {', '.join(data.columns)}")
        return
    
    # Filter for target chromosome if specified
    if target_chromosome:
        data = data[data['CHR'] == target_chromosome]
        print(f"Filtered to {len(data)} variants on chromosome {target_chromosome}")
    
    # Calculate Bonferroni threshold
    bonferroni_threshold = -np.log10(pvalue_threshold / len(data))
    print(f"Bonferroni significance threshold: {bonferroni_threshold:.2f} (-log10 p-value)")
    
    # Parse GFF file
    try:
        cds_features = parse_gff(gff_file, target_chromosome)
        print(f"Loaded {len(cds_features)} CDS features from {gff_file}")
    except Exception as e:
        print(f"Error reading GFF file: {e}")
        return
    
    # Identify significant variants
    significant_variants = data[data['LOG10P'] >= bonferroni_threshold].copy()
    print(f"Found {len(significant_variants)} significant variants")
    
    # Find nearby genes for significant variants
    significant_variants['gene_annotation'] = significant_variants.apply(
        lambda row: find_nearby_genes(row['BP'], cds_features, distance), axis=1
    )
    
    # Create Manhattan plot
    plt.figure(figsize=(12, 6), dpi=300)
    
    # Plot all variants
    plt.scatter(data['BP'], data['LOG10P'], s=10, alpha=0.7, color='blue', label='Variants')
    
    # Highlight significant variants
    if len(significant_variants) > 0:
        plt.scatter(significant_variants['BP'], significant_variants['LOG10P'], 
                    s=30, alpha=1.0, color='red', label='Significant variants')
    
    # Add significance threshold line
    plt.axhline(y=bonferroni_threshold, color='red', linestyle='--', 
                label=f'Bonferroni threshold (p={pvalue_threshold})')
    
    # Add gene annotations for significant variants
    for _, variant in significant_variants.iterrows():
        if variant['gene_annotation'] and len(variant['gene_annotation']) > 0:
            nearby_gene = variant['gene_annotation'][0]  # Get the closest gene
            label = f"{nearby_gene['gene_name']} ({nearby_gene['distance']}bp)"
            
            # Position the label above the point
            plt.annotate(label, 
                        (variant['BP'], variant['LOG10P']),
                        textcoords="offset points", 
                        xytext=(0, 10),
                        ha='center',
                        fontsize=8,
                        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    
    # Set plot labels and title
    chromosome_label = f"Chromosome {target_chromosome}" if target_chromosome else "Chromosome"
    plt.xlabel(f'{chromosome_label} Position (bp)')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot with Gene Annotations')
    
    # Add legend
    plt.legend(loc='upper right')
    
    # Add grid
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Format x-axis to show scientific notation for large numbers
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    print(f"Manhattan plot saved to {output_file}")
    
    # Create a table of significant variants with gene annotations
    if len(significant_variants) > 0:
        results_table = []
        
        for _, variant in significant_variants.iterrows():
            variant_info = {
                'CHR': variant['CHR'],
                'Position': variant['BP'],
                'P-value': variant['P'],
                '-log10(P)': variant['LOG10P']
            }
            
            if variant['gene_annotation'] and len(variant['gene_annotation']) > 0:
                nearby_gene = variant['gene_annotation'][0]
                variant_info.update({
                    'Nearby_gene': nearby_gene['gene_name'],
                    'Gene_product': nearby_gene['product'],
                    'Distance_bp': nearby_gene['distance'],
                    'Strand': nearby_gene['strand']
                })
            else:
                variant_info.update({
                    'Nearby_gene': 'None',
                    'Gene_product': 'None',
                    'Distance_bp': 'NA',
                    'Strand': 'NA'
                })
            
            results_table.append(variant_info)
        
        # Save table to file
        results_df = pd.DataFrame(results_table)
        table_output = output_file.replace('.png', '_significant_variants.tsv')
        results_df.to_csv(table_output, sep='\t', index=False)
        print(f"Table of significant variants saved to {table_output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a Manhattan plot with gene annotations')
    
    parser.add_argument('--manhattan', required=True, 
                        help='Path to the Manhattan plot data file (tab-separated with columns: CHR, BP, P, LOG10P)')
    
    parser.add_argument('--gff', required=True, 
                        help='Path to the GFF annotation file')
    
    parser.add_argument('--chromosome', default=None, 
                        help='Target chromosome to plot (default: plot all chromosomes)')
    
    parser.add_argument('--distance', type=int, default=1000, 
                        help='Maximum distance (bp) to consider a gene as nearby a variant (default: 1000)')
    
    parser.add_argument('--pvalue', type=float, default=0.05, 
                        help='P-value threshold for significance before Bonferroni correction (default: 0.05)')
    
    parser.add_argument('--output', default='manhattan_plot.png', 
                        help='Output file name (default: manhattan_plot.png)')
    
    args = parser.parse_args()
    
    create_manhattan_plot(
        args.manhattan, 
        args.gff, 
        args.output,
        args.chromosome,
        args.distance,
        args.pvalue
    )

#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# Set the font size and style for better readability
plt.style.use('seaborn-v0_8-whitegrid')

# Read the data
try:
    df = pd.read_csv('mortality_manhattan.txt', sep='\t')
except FileNotFoundError:
    print("Error: File 'mortality_manhattan.txt' not found.")
    exit(1)
except pd.errors.ParserError:
    print("Error: File 'mortality_manhattan.txt' is not properly formatted.")
    exit(1)

# Check for empty DataFrame or missing columns
if df.empty:
    print("Error: The input file is empty.")
    exit(1)

required_columns = ['BP', 'LOG10P']
if not all(col in df.columns for col in required_columns):
    print(f"Error: The input file must contain the columns: {', '.join(required_columns)}")
    print(f"Found columns: {', '.join(df.columns)}")
    exit(1)

# Clean the data
df = df.dropna(subset=['LOG10P'])
df = df[np.isfinite(df['LOG10P'])]

# Define colors for alternating chromosomes
colors = ['#1f77b4', '#ff7f0e']  # blue and orange

# Calculate Bonferroni correction threshold
alpha = 0.05
n_tests = len(df)
if n_tests == 0:
    print("Error: No data points found in the input file.")
    exit(1)

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
else:
    print("No significant SNPs found.")

# Add labels and title
plt.xlabel('Position (bp)')
plt.ylabel('-log10(p-value)')
plt.title('Manhattan Plot for Mortality GWAS')

# Set y-axis to start at 0
y_max = max(df['LOG10P'].max() * 1.05, bonferroni_threshold * 1.1)
plt.ylim([0, y_max])

# Add legend
plt.legend(loc='upper right')

# Tighten layout
plt.tight_layout()

# Check if output file exists
if os.path.exists('mortality_manhattan_plot.png'):
    print("Warning: The file 'mortality_manhattan_plot.png' already exists and will be overwritten.")

# Save figure
plt.savefig('mortality_manhattan_plot.png', dpi=600)
print("Plot saved as mortality_manhattan_plot.png")
