from pastml.acr import pastml_pipeline

data = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage_metadata.csv'
tree = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold_50_highburnin.finaltree.nwk'
columns = ['Year', 'Lineage', 'Location']
html = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage1_metadata_tree.html'
html_output = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage1_metadata.html'
pastml_pipeline(data=data, data_sep=',', tree=tree, columns=columns, html_compressed=html_output, html=html, name_column='Year', verbose=True)
