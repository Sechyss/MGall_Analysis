#!/bin/bash

generate_ska_alignment.py --reference ../CheckMbins/R.fna --input List_files_HF_Rlow.txt --out Rlow_gubbins.aln
run_gubbins.py --prefix VA94_all_gubbins_run --date ~/OneDrive/Data/Cambridge_Project/Edited_all_samples_dates.txt --threads 10 --tree-builder iqtree --bootstrap 1000 --seed 1234 --best-model --model-fitter iqtree --recon-with-dates --filter-percentage 50 --extensive-search HF_V94_96HFgenomes.fasta
mask_gubbins_aln.py --aln Rlow_gubbins.aln --gff Rlow_gubbins_run.recombination_predictions.gff --out Rlow_gubbins_run.masked.aln