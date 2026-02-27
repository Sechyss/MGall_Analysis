# BacDating.R — Bayesian molecular-clock dating of the M. gallisepticum phylogeny.
#
# Uses the BactDating package to fit a relaxed-clock Bayesian model to the
# core-genome ML tree, estimating substitution rates and divergence times.
# The root-to-tip regression is visualised prior to Bayesian inference.
#
# Requirements: BactDating, ape
# Input:  core-genome_HF.treefile  (IQ-TREE output; update setwd/path)
# Output: MCMC posterior samples and time-calibrated tree plots

library(BactDating)
library(ape)
setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Pangenome_results_HF/Iqtree_phylogenies/")

# Load the tree
tree=read.tree("core-genome_HF.treefile")
plot(tree)
ape::axisPhylo(backward=F, root.time = 1990)
obsphy <- simobsphy(tree)

# Plot the root-to-tip regression
res <- bactdate(obsphy, 1990:2015)
plot(res, 'treeCI')