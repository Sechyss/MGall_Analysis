#!/bin/bash
# Define the source directory
source_dir="$1"


for file in "$source_dir"/*fna; do
    filename=$(basename "$file")
    newname="${filename%%.fna*}"
    output_dir="/home/albertotr/OneDrive/Data/Cambridge_Project/CheckMbins/${newname}_vibrant"
    # Check if the output directory already exists
  if [ -d "$output_dir" ]; then
    echo "Output directory $output_dir already exists. Skipping..."
    continue
  fi
    python3 ~/Software/Vibrant_db/VIBRANT/VIBRANT_run.py -i /home/albertotr/OneDrive/Data/Cambridge_Project/"$file" -t 10 -d ~/Software/Vibrant_db/VIBRANT/databases -m ~/Software/Vibrant_db/VIBRANT/files -folder "$output_dir"
done