from pastml.acr import pastml_pipeline

data = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage_metadata.csv'
tree = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/VA94_consensus_all_60thres.finaltree.nwk'
columns = ['Year', 'Lineage']
html = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage1_metadata_tree.html'
html_output = '/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/lineage1_metadata.html'
pastml_pipeline(data=data, data_sep=',', tree=tree, columns=columns, html_compressed=html_output, html=html, name_column='Year', verbose=True)
