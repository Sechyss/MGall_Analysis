library(panstripe)
library(ape)
library(patchwork)
set.seed(1234)


setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Pangenome_results_HF/")

# Load the tree and gene presence absence matrix
pa <- read_rtab("gene_presence_absence_filt_pseudo_length_frag.Rtab")
tree <- read.tree("Iqtree_phylogenies/core-genome_HF.contree")

# Run panstripe pipeline
fit <- panstripe(pa, tree, nboot = 100, conf = 0.95)
# Save plots as PNG in the images folder
png("images/pangenome_params.png", res = 600)
plot_pangenome_params(fit)
dev.off()

png("images/pangenome_cumulative.png", res = 600)
plot_pangenome_cumulative(fit)
dev.off()

png("images/pangenome_curve.png", res = 600)
plot_pangenome_curve(fit)
dev.off()


# Tree with presence absence

# look at only those genes that vary
variable_genes <- colnames(pa)[apply(pa, 2, sd) > 0]

plot_tree_pa(tree = tree, pa = pa, genes = variable_genes, label_genes = FALSE, cols = "black")
plot_gain_loss(fit)