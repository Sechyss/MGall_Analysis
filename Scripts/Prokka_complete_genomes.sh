#!/bin/bash

# Define the source directory
source_dir="$1"

# Loop through each subdirectory in the source directory
for dir in "$source_dir"/; do
    folder_name=$(basename "$dir")
    # Find the .fna file in the current folder
    fna_file=$(find "$dir" -maxdepth 1 -name "G*.fna")
     # Check if .gff file exists
    if [ ! -z "$fna_file" ]; then
      prokka --force --outdir "$folder_name"_prokka --prefix "$folder_name" --kingdom Bacteria --gcode 4 --locustag "$folder_name" --genus 'Mycoplasma' --species 'gallisepticum' --usegenus --rnammer "$fna_file";
    fi
done

