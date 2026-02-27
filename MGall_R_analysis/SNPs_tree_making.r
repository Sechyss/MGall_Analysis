# SNPs_tree_making.r — Visualise the phylogenetic tree with ggtree.
#
# Loads a phylogenetic tree and renders it using ggtree with lineage-based
# colouring and metadata annotation for inclusion in the manuscript.
#
# Requirements: ggtree, TDbook
# Input:  Phylogenetic tree file (update path below)
# Output: Annotated phylogenetic tree plot

library(ggtree)
library(TDbook)