#!/bin/bash
#
# SH_comparison.sh — Test tree topology hypotheses using IQ-TREE SH/AU tests.
#
# Concatenates four alignment variants into a single reference alignment,
# builds individual ML trees for each variant, then uses IQ-TREE to
# evaluate the SH (Shimodaira-Hasegawa) and AU (approximately unbiased)
# topology tests across all candidate trees.
#
# Alignment variants compared:
#   60% occupancy threshold  (VA94_consensus_all_trimmed_60threshold.fasta)
#   80% occupancy threshold  (VA94_consensus_all_trimmed_80threshold.fasta)
#   Parsimony-informative    (VA94_consensus_all_trimmed_informative_sites.fasta)
#   No-gaps                  (VA94_consensus_all_trimmed_nogaps.fasta)
#
# Usage:
#   bash SH_comparison.sh
#   (Run from the directory containing all four alignment FASTA files)
#
# Requirements: seqkit, iqtree (>=2.0), conda (for environment switching)


seqkit concat -o VA94_concatenated_forcomparison.fasta VA94_consensus_all_trimmed_60threshold.fasta VA94_consensus_all_trimmed_80threshold.fasta VA94_consensus_all_trimmed_informative_sites.fasta VA94_consensus_all_trimmed_nogaps.fasta 
conda activate Gubbins_env
iqtree -s ../VA94_consensus_all_trimmed_60threshold.fasta -m GTR+F+I+G4 -T AUTO -B 1000 --prefix VA94_60thr && iqtree -s ../VA94_consensus_all_trimmed_80threshold.fasta -m GTR+F+I+G4 -T AUTO -B 1000 --prefix VA94_80thr && iqtree -s ../VA94_consensus_all_trimmed_informative_sites.fasta -m GTR+F+I+G4 -T AUTO -B 1000 --prefix VA94_infopars && iqtree -s ../VA94_consensus_all_trimmed_nogaps.fasta -m GTR+F+I+G4 -T AUTO -B 1000 --prefix VA94_nogaps

#Combine all trees into a single file
cat VA94_60thr.treefile VA94_80thr.treefile VA94_infopars.treefile VA94_nogaps.treefile > Combined_trees.txt
iqtree -s ../VA94_concatenated_forcomparison.fasta -z Combined_trees.txt -T AUTO -m TEST -zb 10000 -zw -au -n 0 -safe -redo --prefix SH_alltrees_VA94_all