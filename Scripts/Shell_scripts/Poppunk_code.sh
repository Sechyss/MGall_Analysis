#!/bin/bash

# Run PopPUNK on the VA94 dataset
# shellcheck disable=SC2164
cd /mnt/c/Users/at991/OneDrive\ -\ University\ of\ Exeter/Data/Cambridge_Project/PopPUNK/
poppunk --create-db --output mgall_60threshold_VA94 --r-files List_files_VA94-Trimmed.txt --threads 8
poppunk --fit-model lineage --ref-db mgall_60threshold_VA94 --ranks 1,2,3,4,5 --overwrite