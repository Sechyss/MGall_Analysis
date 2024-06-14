#!/bin/bash

# shellcheck disable=SC2164
cd /home/albertotr/OneDrive/Data/Cambridge_Project/SPADES_denovo_SRA/

panaroo -i SR*_prokka/*.gff -o ../pangenome_results/ --clean-mode strict -a core --aligner clustal --core_threshold 0.98 -t 10