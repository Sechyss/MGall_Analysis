library(BactDating)
library(ape)
setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Pangenome_results_HF/Iqtree_phylogenies/")

# Load the tree
tree=simcoaltree(1990:2010)
plot(tree)
ape::axisPhylo(backward=F)
