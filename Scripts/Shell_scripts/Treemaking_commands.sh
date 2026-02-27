#!/bin/bash
#
# Treemaking_commands.sh — Build phylogenetic trees with IQ-TREE, FastTree, and RAxML.
#
# Reference commands for constructing maximum-likelihood trees from the
# Rlow SNP alignment using three different tree-building tools.
# Run each command individually or together depending on requirements.
#
# Tools and key settings:
#   IQ-TREE   Automatic model selection (-m TEST), 1000 ultrafast bootstraps,
#             time-calibration using sample dates (--date)
#   FastTree  GTR+CAT model, 1000 bootstrap replicates
#   RAxML     GTRGAMMA model, rapid bootstrap (100 replicates)
#
# Usage:
#   bash Treemaking_commands.sh
#   (Run from the directory containing the masked alignment FASTA)
#
# Requirements: iqtree (>=2.0), FastTree, raxmlHPC-PTHREADS

iqtree --prefix Rlowsnps_stripped_nogaps -s Rlow_consensus_snps_trimmed_nogaps.masked.fasta -m TEST -T AUTO -B 1000 --date ~/OneDrive/Data/Cambridge_Project/Edited_Lucy_samples_dates.txt

fasttree -nt -gtr -boot 1000 < Rlow_consensus_snps_trimmed_nogaps.masked.fasta > Fasttree_Rlow_snps_Lucy_trimmed.newick

raxmlHPC-PTHREADS -T 10 -s Rlow_consensus_snps_trimmed_nogaps.masked.fasta -f a -n RaxML_Rlowsnps_stripped_nogaps -m GTRGAMMA4 -c -F -k -b 1234 --bootstop-perms=100