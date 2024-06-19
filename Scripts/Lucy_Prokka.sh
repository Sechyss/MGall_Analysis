#!/bin/bash
# Define the source directory
source_dir="$1"

# List all files in the directory and process each one
for dir in "$source_dir"/*_assembly/; do
  # Extract the filename from the full path
  dirname=$(basename "$dir")

  # Cut the filename after the string 'sickle'
  newname="${dirname%%_nophi_assembly*}"

  prokka --force --outdir "$source_dir"/"$newname"_prokka --prefix "$newname" --kingdom Bacteria --gcode 4 --locustag "$newname" --genus 'Mycoplasma' --species 'gallisepticum' --usegenus --rnammer "$dir"contigs.fasta; done