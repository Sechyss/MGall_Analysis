#!/bin/bash
#
# Panaroo_pipeline.sh — Build the M. gallisepticum pangenome with Panaroo.
#
# Runs Panaroo on all house-finch GFF annotation files to construct a
# high-quality pangenome graph and gene presence/absence matrix.
#
# Key parameters:
#   --clean-mode strict       Aggressively remove likely contaminants / artefacts
#   --codon-table 4           Mycoplasma/Spiroplasma genetic code (UGA = Trp)
#   --remove-invalid-genes    Discard genes with internal stop codons (table 4)
#   -a core                   Align core-genome genes
#   --core_threshold 0.98     Gene present in ≥98% of isolates = core gene
#   --aligner mafft           Use MAFFT for multi-sequence alignment
#   --refind_prop_match 0.5   Minimum proportion for re-finding missing genes
#   --search_radius 1000      Search radius (bp) for missing gene re-finding
#
# Usage:
#   bash Panaroo_pipeline.sh
#   (Update input GFF path and output directory before running)
#
# Requirements: panaroo (>=1.3), mafft

 /mnt/c/Users/at991/OneDrive\ -\ University\ of\ Exeter/Data/Cambridge_Project/GFFs_genomes/HF_gffs/*.gff -o Pangenome_results_HF/ --clean-mode strict --codon-table 4 --remove-invalid-genes -a core --core_threshold 0.98 -t 10 --aligner mafft --refind_prop_match 0.5 --search_radius 1000