#!/bin/bash

cd ~/OneDrive/Data/Ipoutcha_motility/BLAST_db/ || exit

for file in ~/OneDrive/Data/Ipoutcha_motility/BLAST_db/*fna; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  if [ -e "$filename".aln ]; then
        echo "Output file $filename already exists. Skipping $filename."
        continue
  fi
  muscle -align "$file" -output "$filename".aln
done

conda init
conda activate Gubbins_env
cd ./Phylogenetic_trees/ || exit

for file in ~/OneDrive/Data/Ipoutcha_motility/BLAST_db/*aln; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  iqtree --prefix Candidate_"$filename" -s ../"$filename".aln -b 1000 -m TEST -T AUTO
done