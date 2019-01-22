import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import plotly.plotly as py
import plotly.graph_objs as go
from plotly.offline.offline import matplotlib

# plot dimensions
matplotlib.rcParams['figure.figsize'] = [15, 20] # for square canvas
matplotlib.rcParams['figure.subplot.left'] = 0
matplotlib.rcParams['figure.subplot.bottom'] = 0
matplotlib.rcParams['figure.subplot.right'] = 1
matplotlib.rcParams['figure.subplot.top'] = 1

# read flagged interactions: colums 0 and 1 as nodes, rows as interactions
data = pd.read_csv('../outputs/flagged_interactions.tsv', sep="\t", header=None)
interactions = [tuple([row[0], row[1]]) for row in data.values.tolist() if row[2] == 1]
drug_nodes = list(set([row[0] for row in interactions]))
gene_nodes = list(set([row[1] for row in interactions]))

# create bipartite graph
G = nx.Graph()
# Add nodes with the node attribute "bipartite"
G.add_nodes_from(drug_nodes, bipartite=0)
G.add_nodes_from(gene_nodes, bipartite=1)
# Add edges only between nodes of opposite node sets
G.add_edges_from(interactions)
nx.is_connected(G)

# distribute node positions according to their bipartite 'class' 0-drugs, 1-genes
X, Y = bipartite.sets(G)
pos = dict()
pos.update((n, (0, i*10)) for i, n in enumerate(X))
pos.update((n, (0.5, i*100)) for i, n in enumerate(Y))
num_edges = G.number_of_edges()
num_nodes = G.number_of_nodes()

# draw in the positions
nx.draw(G, pos=pos, with_labels=True,
        edge_color="Green",
        edge_cmap=plt.get_cmap('Greens'),
        node_color="Red",
        cmap=plt.get_cmap('Reds'))

plt.savefig("../outputs/NNETGraph/bipartite.png")
