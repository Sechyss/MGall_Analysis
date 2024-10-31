import allel
import numpy as np

vcf_file = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_Rlow/Only_SNPs/A1.bam.vcf.gz'
genotypes = allel.read_vcf(vcf_file)['calldata/GT']

# Identify segregating sites (polymorphic sites)
# A site is segregating if it has more than one allele present across samples
segregating_sites = (genotypes.max(axis=1) > 0) & (genotypes.min(axis=1) == 0)

# Print segregating site positions (True indicates segregating site)
for row in segregating_sites:
    if True in row:
        print(row)