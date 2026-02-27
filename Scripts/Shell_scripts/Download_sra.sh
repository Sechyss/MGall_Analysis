#!/bin/bash
#
# Download_sra.sh — Batch download paired-end FASTQ files from EBI FTP.
#
# Downloads R1 and R2 FASTQ.gz files for three sets of SRA accessions
# (SRR628xxxx house-finch isolates and two SRR131xxxx isolates) from the
# EBI FTP mirror using wget with -nc (no-clobber) to skip existing files.
# Downloads run in parallel background processes within each batch.
#
# Usage:
#   bash Download_sra.sh
#
# Requirements: wget
# Output directory: current working directory

 "SRR6289955" "SRR6289956" "SRR6289957" "SRR6289958" "SRR6289959" "SRR6289960" "SRR6289961" "SRR6289962" "SRR6289963" "SRR6289964" "SRR6289965" "SRR6289966" "SRR6289970" "SRR6289971" "SRR6289973" "SRR6289977" "SRR6289978" "SRR6289979" "SRR6289980" "SRR6289982" "SRR6289983" "SRR6289984" "SRR6289985" "SRR6289986" "SRR6289987" "SRR6289988" "SRR6289989" "SRR6289990" "SRR6289991" "SRR6289992" "SRR6289993" "SRR6289994" "SRR6289995" "SRR6289996" "SRR6289997" "SRR6289998" "SRR6289999")
runs2=("SRR6290000" "SRR6290001" "SRR6290002" "SRR6290003" "SRR6290004" "SRR6290005" "SRR6290006" "SRR6290007" "SRR6290008")
runs3=("SRR13148778" "SRR13148779")

extract_sra(){
  id_run=$1
  batch=$2
  echo "Downloading ${id_run}"
  
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$batch"/008/"$id_run"/"$id_run"_1.fastq.gz
  wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"$batch"/008/"$id_run"/"$id_run"_2.fastq.gz
}

export -f extract_sra
for run_id in "${runs2[@]}"; do
  (extract_sra "$run_id" "SRR629") & done

wait


