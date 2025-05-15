#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import argparse
from matplotlib.ticker import ScalarFormatter
import os
from BCBio import GFF
from Bio import SeqIO


def parse_gff(gff_file, target_chromosome=None):
    """
    Parse a GFF file and extract CDS features with gene information using BCBio.
    
    Args:
        gff_file (str): Path to the GFF file
        target_chromosome (str, optional): If provided, only features from this chromosome are extracted
        
    Returns:
        list: List of dictionaries containing CDS features with their positions and annotations
    """
    try:
        cds_features = []
        with open(gff_file, 'r') as f:
            for record in GFF.parse(f):
                if target_chromosome and record.id != target_chromosome:
                    continue
                for feature in record.features:
                    if feature.type == 'CDS':
                        gene_name = feature.qualifiers.get('gene', ['Unknown'])[0]
                        product = feature.qualifiers.get('product', ['Unknown'])[0]
                        strand = '+' if feature.location.strand == 1 else '-'  # Updated to use .location.strand
                        cds_features.append({
                            'chromosome': record.id,
                            'start': int(feature.location.start) + 1,  # Convert to 1-based indexing
                            'end': int(feature.location.end),
                            'strand': strand,
                            'gene_name': gene_name,
                            'product': product
                        })
        # Debugging: Print the number of CDS features parsed
        print(f"Parsed {len(cds_features)} CDS features from GFF file.")
        return cds_features
    except FileNotFoundError:
        print(f"Error: GFF file '{gff_file}' not found.")
        return []
    except Exception as e:
        print(f"Error parsing GFF file: {e}")
        return []


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
    
    # Debugging: Print nearby genes for the given position
    if not nearby_genes:
        print(f"No nearby genes found for position {position}.")
    else:
        print(f"Nearby genes for position {position}: {nearby_genes}")
    
    return nearby_genes if nearby_genes else [{'gene_name': 'Unknown', 'product': 'None', 'distance': 0, 'strand': 'NA'}]


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
    try:
        data = pd.read_csv(manhattan_data, sep='\t')
        print(f"Loaded {len(data)} variants from {manhattan_data}")
    except FileNotFoundError:
        print(f"Error: Manhattan data file '{manhattan_data}' not found.")
        return
    except Exception as e:
        print(f"Error reading Manhattan data file: {e}")
        return

    if data.empty:
        print("Error: The Manhattan data file is empty.")
        return

    required_columns = ['CHR', 'BP', 'P', 'LOG10P']
    if not all(col in data.columns for col in required_columns):
        print(f"Error: Manhattan data file must contain columns: {', '.join(required_columns)}")
        print(f"Found columns: {', '.join(data.columns)}")
        return

    data = data.dropna(subset=['LOG10P'])
    data = data[np.isfinite(data['LOG10P'])]

    if len(data) == 0:
        print("Error: No data points found in the Manhattan data file.")
        return

    bonferroni_threshold = -np.log10(pvalue_threshold / len(data))
    print(f"Bonferroni significance threshold: {bonferroni_threshold:.2f} (-log10 p-value)")

    cds_features = parse_gff(gff_file, target_chromosome)
    if not cds_features:
        return

    significant_variants = data[data['LOG10P'] >= bonferroni_threshold].copy()
    if significant_variants.empty:
        print("No significant variants found.")
        return

    significant_variants['gene_annotation'] = significant_variants.apply(
        lambda row: find_nearby_genes(row['BP'], cds_features, distance), axis=1
    )

    plt.figure(figsize=(12, 6), dpi=300)
    plt.scatter(data['BP'], data['LOG10P'], s=10, alpha=0.7, color='blue', label='Variants')
    plt.axhline(y=bonferroni_threshold, color='red', linestyle='--', 
                label=f'Bonferroni threshold (p={pvalue_threshold})')

    if len(significant_variants) > 0:
        plt.scatter(significant_variants['BP'], significant_variants['LOG10P'], 
                    s=30, alpha=1.0, color='red', label='Significant variants')

    for _, variant in significant_variants.iterrows():
        if variant['gene_annotation'] and len(variant['gene_annotation']) > 0:
            nearby_gene = variant['gene_annotation'][0]
            label = f"{nearby_gene['gene_name']} ({nearby_gene['distance']}bp)"
            plt.annotate(label, 
                         (variant['BP'], variant['LOG10P']),
                         textcoords="offset points", 
                         xytext=(0, 10),
                         ha='center',
                         fontsize=8,
                         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))

    plt.xlabel('Position (bp)')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot with Gene Annotations')
    plt.legend(loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    if os.path.exists(output_file):
        print(f"Warning: The file '{output_file}' already exists and will be overwritten.")
    plt.savefig(output_file, dpi=600)
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
