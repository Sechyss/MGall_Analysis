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