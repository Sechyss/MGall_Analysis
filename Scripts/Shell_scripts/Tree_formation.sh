#!/bin/bash
#
# Tree_formation.sh — Build per-gene phylogenetic trees for motility candidates.
#
# For each motility candidate .fna FASTA file in the BLAST database directory,
# this script aligns sequences with MUSCLE and then builds a maximum-likelihood
# tree with IQ-TREE (automatic model selection, 1000 bootstraps).
# Already-aligned files are skipped.
#
# Usage:
#   bash Tree_formation.sh
#   (Update the base directory path before running)
#
# Requirements: muscle (>=5), iqtree (>=2.0), conda (for environment activation)
# Outputs:
#   <gene>.aln           MUSCLE multiple sequence alignment
#   Phylogenetic_trees/  IQ-TREE output files per gene



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