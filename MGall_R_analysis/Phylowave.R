# ---- Loading basic data ----

source(file = '2_Functions/2_1_Index_computation_20240909.R')
source(file = '2_Functions/2_2_Lineage_detection_20240909.R')
source(file = '2_Functions/2_3_Lineage_fitness_20240909.R')

library(ape, quiet = T); library(phytools, quiet = T); library(stringr, quiet = T)
library(MetBrewer, quiet = T); library(parallel, quiet = T); library(mgcv, quiet = T)
library(cowplot, quiet = T); library(ggplot2, quiet = T); library(ggtree, quiet = T);
library(cmdstanr, quiet = T); library(binom, quiet = T)

setwd("C:/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs")

timed_tree = read.nexus('./Rlowsnps_stripped_nogaps.timetree.nex')

genome_length = 996422


#---- Names all sequences ----
names_seqs = timed_tree$tip.label
n_seq = length(names_seqs)
## Collection times of all sequences
times_seqs = as.numeric(sapply(names_seqs, function(x)tail(str_split(x, pattern = '/')[[1]],2)[1]))
## Nextstrain clades of all sequences
clades_seqs = sapply(names_seqs, function(x)tail(str_split(x, pattern = '/')[[1]],1))

## Pairwise distance matrix
genetic_distance_mat = dist.nodes.with.names(timed_tree)

## Time for each internal node

nroot = length(timed_tree$tip.label) + 1 ## Root number
distance_to_root = genetic_distance_mat[nroot,]
root_height = times_seqs[which(names_seqs == names(distance_to_root[1]))] - distance_to_root[1]
nodes_height = root_height + distance_to_root[n_seq+(1:(n_seq-1))]

