#!/bin/bash

panaroo -i /mnt/c/Users/at991/OneDrive\ -\ University\ of\ Exeter/Data/Cambridge_Project/GFFs_genomes/HF_gffs/*.gff -o Pangenome_results_HF/ --clean-mode strict --codon-table 4 --remove-invalid-genes -a core --core_threshold 0.98 -t 10 --aligner mafft --refind_prop_match 0.5 --search_radius 1000