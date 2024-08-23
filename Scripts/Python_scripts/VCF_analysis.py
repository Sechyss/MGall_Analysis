import pysam

vcf_reader = pysam.VariantFile('/home/albertotr/OneDrive/Data/'
                               'Cambridge_Project/Mapped_output_VA94_7994_1_7P/A1.bam.vcf.gz', 'r')

for record in vcf_reader:
