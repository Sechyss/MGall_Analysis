#!/bin/bash

cd ~/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/VCF_files || exit

for file in ~/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BAMFiles/*bam; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  if [ -e "$filename".fasta ]; then
        echo "Output file $filename already exists. Skipping $filename."
        continue
  fi
  bcftools mpileup -f ~/OneDrive/Data/Cambridge_Project/CheckMbins/VA94_7994_1_7P.fna --threads 20 "$file" | bcftools call --threads 20 --ploidy 1 -m -Ob -o calls.bcf
  bcftools view --threads 20 -V indels -i 'QUAL>=10' calls.bcf > "$filename".vcf
  bgzip -f "$filename".vcf
  bcftools index --threads 20 "$filename".vcf.gz
  bcftools consensus -a - --mark-snv lc -f ~/OneDrive/Data/Cambridge_Project/CheckMbins/VA94_7994_1_7P.fna "$filename".vcf.gz -o "$filename".fasta
done
