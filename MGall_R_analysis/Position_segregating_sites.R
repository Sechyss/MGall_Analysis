# Load the ape package
library(ape)

setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/")

# Path to the DNA alignment file (FASTA format or other readable by ape)
alignment_file <- "Rlow_consensus_strains_onlyLucy.fasta"
lucy_alignment <- "../../VCF_trees/Myco_py.fasta"
# Read the DNA alignment
alignment <- read.dna(alignment_file, format = "fasta")
alignment_lucy <- read.dna(lucy_alignment, format = "fasta")
# Identify segregating sites
segregating_sites <- seg.sites(alignment)
segregating_sites_lucy <- seg.sites(alignment_lucy)

# Convert lists to vectors
vec1 <- unlist(segregating_sites)
vec2 <- unlist(segregating_sites_lucy)

# Find common elements
common_elements <- intersect(vec1, vec2)
print(common_elements)

# Find elements unique to each list
unique_to_list1 <- setdiff(vec1, vec2)
unique_to_list2 <- setdiff(vec2, vec1)
