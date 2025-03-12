from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from collections import Counter
import io

# Load data
snps_lucy = pd.read_csv('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_VA94_7994_1_7P/All_mutations_matrix.csv', index_col=0)
snps_sra = pd.read_excel('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/All_mutation_matrix.xlsx', index_col=0, sheet_name='Sheet1', engine='openpyxl')

lucy_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/Lucy_replacements.pickle', 'rb'))
sra_replacement = pickle.load(open('/home/albertotr/OneDrive/Data/Cambridge_Project/SRA_replacements.pickle', 'rb'))

snps_lucy['Sample'] = snps_lucy['Sample'].replace(lucy_replacement)
snps_sra['Sample'] = snps_sra['Sample'].replace(sra_replacement)

all_snps = pd.concat([snps_lucy, snps_sra], axis=0)

non_synonymous = all_snps[all_snps['Effect'].str.contains('STOP') |(all_snps['Effect'].str.contains('NON_SYNONYMOUS')) | (all_snps['Effect'].str.contains('LOST'))]

# Load tree
tree_file = Tree('/home/albertotr/OneDrive/Data/Cambridge_Project/Mapped_output_SRA_VA94/BEAST/Final_Run/VA94_consensus_all_trimmed_60threshold.finaltree.nwk')

# Function to create pie chart
def create_pie_chart(mutations):
    counts = Counter(mutations)
    labels, sizes = zip(*counts.items())
    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels, autopct='%1.1f%%')
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    plt.close(fig)
    return buf

# Add pie charts to nodes
for node in tree_file.traverse():
    if not node.is_leaf():
        # Get mutations for this node
        node_mutations = non_synonymous[non_synonymous['Node'] == node.name]['Mutation']
        if not node_mutations.empty:
            pie_chart = create_pie_chart(node_mutations)
            pie_face = faces.ImgFace(pie_chart, width=100, height=100)
            node.add_face(pie_face, column=0, position="branch-right")

# Render tree
ts = TreeStyle()
ts.show_leaf_name = True
tree_file.show(tree_style=ts)


