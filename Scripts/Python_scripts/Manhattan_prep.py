#!/usr/bin/env python3
import math
import os

os.chdir('/home/albertotr/OneDrive/Data/Cambridge_Project/GWAS/')

# Input and output file paths
input_file = "swelling_SNPs.txt"
output_file = "swelling_manhattan.txt"

# Function to extract chromosome and position from variant ID
def parse_variant(variant):
    # Format: chromosome_position_ref_alt
    parts = variant.split('_')
    chromosome = parts[0]
    position = int(parts[1]) if parts[1].isdigit() else parts[1]
    return chromosome, position

# Process the file
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Write header to output file
    outfile.write("CHR\tBP\tP\tLOG10P\n")
    
    # Skip header line from input
    header = infile.readline()
    
    # Process each line
    for line in infile:
        columns = line.strip().split('\t')
        
        # Extract variant info
        variant = columns[0]
        try:
            chromosome, position = parse_variant(variant)
            
            # Get p-value (3rd column, index 2)
            p_value_str = columns[2]
            p_value = float(p_value_str)
            
            # Calculate -log10(p)
            if p_value > 0:  # Avoid log(0)
                log10p = -math.log10(p_value)
            else:
                log10p = 0  # or some maximum value
                
            # Write to output file
            outfile.write(f"{chromosome}\t{position}\t{p_value}\t{log10p}\n")
            
        except Exception as e:
            print(f"Error processing line: {line.strip()}")
            print(f"Error details: {e}")

print(f"Processing complete. Output written to {output_file}")