# Load the ape package
library(ape)

setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/")

# Path to the DNA alignment file (FASTA format or other readable by ape)
alignment_file <- "Rlow_consensus_strains_onlyLucy.fasta"
# Read the DNA alignment
alignment <- read.dna(alignment_file, format = "fasta")

# Identify segregating sites
segregating_sites <- seg.sites(alignment)

# Report positions in the alignment that contain segregating sites
if (length(segregating_sites) > 0) {
  cat("Positions with segregating sites:\n")
  print(segregating_sites)
} else {
  cat("No segregating sites found.\n")
}
