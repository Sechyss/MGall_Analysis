#!/bin/bash
#
# Gubbins_pipeline.sh — Detect and mask recombination in a whole-genome alignment.
#
# Runs Gubbins on the combined house-finch alignment to identify recombinant
# regions and mask them prior to phylogenetic analysis.  After recombination
# detection the masked alignment is produced with mask_gubbins_aln.py.
#
# Key options:
#   --filter-percentage 50   Remove sites present in >50% of recombination blocks
#   --extensive-search       More thorough recombination detection
#   --recon-with-dates       Time-aware ancestral reconstruction
#   --tree-builder iqtree    Use IQ-TREE for internal tree steps
#   --bootstrap 1000         Bootstrap replicates for internal trees
#
# Usage:
#   bash Gubbins_pipeline.sh
#   (Run from the directory containing HF_V94_96HFgenomes.fasta)
#
# Requirements: gubbins (run_gubbins.py, mask_gubbins_aln.py), IQ-TREE

run_gubbins.py --prefix VA94_all_gubbins_run --date ~/OneDrive/Data/Cambridge_Project/Edited_all_samples_dates.txt --threads 10 --tree-builder iqtree --bootstrap 1000 --seed 1234 --best-model --model-fitter iqtree --recon-with-dates --filter-percentage 50 --extensive-search HF_V94_96HFgenomes.fasta
mask_gubbins_aln.py --aln Rlow_gubbins.aln --gff Rlow_gubbins_run.recombination_predictions.gff --out Rlow_gubbins_run.masked.aln