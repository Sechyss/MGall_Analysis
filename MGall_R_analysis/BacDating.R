library(BactDating)
library(ape)
setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Pangenome_results_HF/Iqtree_phylogenies/")

# Load the tree
tree=read.tree("core-genome_HF.timetree.nex")
plot(tree)
ape::axisPhylo(backward=F, root.time = 1990)
