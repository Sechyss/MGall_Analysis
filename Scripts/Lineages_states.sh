#usr/bin/bash

# This script takes an annotated phylogenetic tree based on lineages and generate the states of each lineage over time of the tree



.\applauncher.bat StateT -t Lineage -i C:\Users\at991\Documents\VA94_consensus_all_trimmed_60threshold_50_highburnin.finaltree.newick -burnin 0 -out C:\Users\at991\Documents\VA94_consensus_all_60thres.finaltreeLineages_statetransition.txt -resolution 59

# -resolution 59 is the number of time points in time of the tree, to generate the desired length of the table
# -burnin 0 is the number of time points to be ignored in the tree
# -t Lineage is the annotation variable of the tree
# -i is the input tree
# -out is the output file
# The output file is a table with the states of each lineage over time of the tree