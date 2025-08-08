#!/bin/bash

runs=("SRR6289954" "SRR6289955" "SRR6289956" "SRR6289957" "SRR6289958" "SRR6289959" "SRR6289960" "SRR6289961" "SRR6289962" "SRR6289963" "SRR6289964" "SRR6289965" "SRR6289966" "SRR6289970" "SRR6289971" "SRR6289973" "SRR6289977" "SRR6289978" "SRR6289979" "SRR6289980" "SRR6289982" "SRR6289983" "SRR6289984" "SRR6289985" "SRR6289986" "SRR6289987" "SRR6289988" "SRR6289989" "SRR6289990" "SRR6289991" "SRR6289992" "SRR6289993" "SRR6289994" "SRR6289995" "SRR6289996" "SRR6289997" "SRR6289998" "SRR6289999" "SRR6290000" "SRR6290001" "SRR6290002" "SRR6290003" "SRR6290004" "SRR6290005" "SRR6290006" "SRR6290007" "SRR6290008")

for accession_id in "${runs[@]}"; do
python3 spades.py -1 /home/albertotr/OneDrive/Data/Cambridge_Project/SRA_reads/"$accession_id"_sickle_R1.fastq -2 /home/albertotr/OneDrive/Data/Cambridge_Project/SRA_reads/"$accession_id"_sickle_R2.fastq -s /home/albertotr/OneDrive/Data/Cambridge_Project/SRA_reads/"$accession_id"_sickle_single.fastq -o /home/albertotr/OneDrive/Data/Cambridge_Project/SPADES_denovo_SRA/"$accession_id"_assembly
done