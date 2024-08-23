#!/bin/bash

cd Mapped_output_Rlow/ || exit

for file in ~/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/*bam; do
  # Extract the filename from the full path
  filename=$(basename "$file")
  if [ -e "$filename".fasta ]; then
        echo "Output file $filename already exists. Skipping $filename."
        continue
  fi
  bcftools mpileup -f ~/OneDrive/Data/Cambridge_Project/CheckMbins/R.fna --threads 20 "$file" | bcftools call --threads 20 --ploidy 1 -mv -Ob -o calls.bcf
  bcftools view --threads 20 -V indels -i 'QUAL>=10' calls.bcf > "$filename".vcf
  #sed -i 's/AE015450\.2/AE015450/g' "$filename".vcf
  bgzip -f "$filename".vcf
  bcftools index --threads 20 "$filename".vcf.gz
  bcftools consensus -f ~/OneDrive/Data/Cambridge_Project/CheckMbins/R.fna "$filename".vcf.gz > "$filename".fasta
done

