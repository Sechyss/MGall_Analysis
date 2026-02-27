#!/bin/bash
#
# BCFtools.sh — Generate consensus FASTA sequences from BAM files.
#
# For each BAM file in the mapped-output directory, this script:
#   1. Runs bcftools mpileup + call to produce a VCF of variant calls.
#   2. Filters to SNPs with QUAL >= 10 (indels excluded).
#   3. Compresses and indexes the VCF with bgzip/bcftools index.
#   4. Applies the variants to the reference with bcftools consensus to
#      generate a per-sample FASTA consensus sequence.
#
# Already-processed samples (existing .fasta output) are skipped.
#
# Usage:
#   bash BCFtools.sh
#
# Requirements: bcftools (>=1.9), bgzip, samtools
# Reference: VA94_7994_1_7P.fna

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
