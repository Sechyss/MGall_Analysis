# Panstripe.R — Pangenome gene gain/loss analysis with Panstripe.
#
# Fits Panstripe models to the M. gallisepticum pangenome presence/absence
# matrix to quantify rates of gene gain and loss along the phylogeny and
# visualises the results as stripe plots with bootstrap confidence intervals.
#
# Requirements: panstripe, ape, patchwork
# Input:  Panaroo presence/absence Rtab and ML phylogeny (update paths below)
# Output: PNG stripe plots (gene gain/loss rates per lineage)

library(panstripe)
library(ape)
library(patchwork)
set.seed(1234)

setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Pangenome_results_HF/")

# Load the tree and gene presence absence matrix
pa <- read_rtab("Edited_gene_presence_absence_filt_pseudo_length_frag_Copy.Rtab")
tree <- read.tree("Iqtree_phylogenies/core-genome_HF.contree")

# Compare the index of table pa and the tree labels
pa_index <- rownames(pa)
tree_labels <- tree$tip.label

# Find mismatches
mismatches <- setdiff(pa_index, tree_labels)
if (length(mismatches) > 0) {
    cat("Mismatches found:\n")
    print(mismatches)
} else {
    cat("No mismatches found.\n")
}

# Run panstripe pipeline
fit <- panstripe(pa, tree, nboot = 100, conf = 0.95)

# Save plots as PNG in the images folder with A4 size dimensions
png("images/pangenome_params.png", width = 4960, height = 7016, res = 600)
plot_pangenome_params(fit)
dev.off()

png("images/pangenome_cumulative.png", width = 4960, height = 7016, res = 600)
plot_pangenome_cumulative(fit)
dev.off()

png("images/pangenome_curve.png", width = 4960, height = 7016, res = 600)
plot_pangenome_curve(fit)
dev.off()

# Tree with presence absence

# look at only those genes that vary
variable_genes <- colnames(pa)[apply(pa, 2, sd) > 0]

png("images/presence_absence_tree.png", width = 7016, height = 7016, res = 600)
plot_tree_pa(tree = tree, pa = pa, genes = variable_genes, label_genes = FALSE, cols = "black")
dev.off()

png("images/plot_gain_loss_tree.png", width = 7016, height = 7016, res = 600)
plot_gain_loss(fit)
dev.off()

png("images/accessory_size.png", width = 4960, height = 7016, res = 600)
plot_acc(pa)
dev.off()
