#!/bin/bash -l

# Define the source directory
source_dir="$1"

for file in "$source_dir"/*fna; do
    filename=$(basename "$file")
    newname="${filename%%.fna*}"
    output_dir="/home/albertotr/OneDrive/Data/Cambridge_Project/Padloc_results/"
    output_file="/home/albertotr/OneDrive/Data/Cambridge_Project/Padloc_results/${newname}_padloc.csv"
    if [ -e "$output_file" ]; then
        echo "Output file $newname already exists. Skipping $newname."
        continue
    fi
    padloc --fna /home/albertotr/OneDrive/Data/Cambridge_Project/"$file" --cpu 10 --outdir "$output_dir"
done

