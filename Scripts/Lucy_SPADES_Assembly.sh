#!/bin/bash
# Define the source directory
source_dir="$1"

# List all files in the directory and process each one
for file in "$source_dir"/*_R1.fastq.gz; do
  # Extract the filename from the full path
  filename=$(basename "$file")

  # Cut the filename after the string 'sickle'
  newname="${filename%%sickle*}"

  # Define the output directory
  output_dir="/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_SPADES_denovo/${newname}assembly"

  # Check if the output directory already exists
  if [ -d "$output_dir" ]; then
    echo "Output directory $output_dir already exists. Skipping..."
    continue
  fi

  python3 spades.py -1 /home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/"$newname"sickle_R1.fastq.gz -2 /home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/"$newname"sickle_R2.fastq.gz -s /home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_reads/"$newname"sickle_single.fastq.gz -o "$output_dir"
done
