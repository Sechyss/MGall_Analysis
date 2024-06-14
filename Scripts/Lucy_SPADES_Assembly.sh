#!/bin/bash
# Define the source directory
source_dir="$1"

# List all files in the directory and process each one
for file in "$source_dir"/*_R1.fastq.gz; do
  # Extract the filename from the full path
  filename=$(basename "$file")

  # Cut the filename after the string 'sickle'
  newname=$(echo "$filename" | sed 's/sickle*//')

  python3 spades.py -1 /home/albertotr/OneDrive/Data/Cambridge_Project/SRA_reads/"$newname"_R1.fastq -2 /home/albertotr/OneDrive/Data/Cambridge_Project/SRA_reads/"$newname"_R2.fastq -s /home/albertotr/OneDrive/Data/Cambridge_Project/SRA_reads/"$newname"_single.fastq -o /home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_SPADES_denovo/"$newname"_assembly
done
