#!/bin/bash

generate_ska_alignment.py --reference ../CheckMbins/R.fna --input List_files_HF_Rlow.txt --out Rlow_gubbins.aln
run_gubbins.py --prefix Rlow_gubbins_run --date ~/OneDrive/Data/Cambridge_Project/Lucy_samples_dates.txt --threads 10 --tree-builder iqtree --bootstrap 1000 --seed 1234 --best-model --model-fitter iqtree --recon-with-dates --filter-percentage 50 --extensive-search Rlow_gubbins.aln
mask_gubbins_aln.py --aln Rlow_gubbins.aln --gff Rlow_gubbins_run.recombination_predictions.gff --out Rlow_gubbins_run.masked.aln