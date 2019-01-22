import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('../outputs/interactions.tsv', sep="\t", header=None)
drug_nodes = list(set(data[0].tolist()))
gene_nodes = list(set(data[1].tolist()))
gene_breast_nodes = list(set(data[2].tolist()))
interaction_edges = data.values.tolist()
interaction_edges = [tuple(l) for l in interaction_edges]

B = nx.Graph()
# Add nodes with the node attribute "bipartite"
B.add_nodes_from(drug_nodes, bipartite=0)
B.add_nodes_from(gene_nodes, bipartite=1)
# Add edges only between nodes of opposite node sets
B.add_edges_from(interaction_edges)
nx.is_connected(B)

bottom_nodes, top_nodes = bipartite.sets(B)

X, Y = bipartite.sets(B)
pos = dict()
pos.update((n, (1, i)) for i, n in enumerate(X))  # put nodes from X at x=1
pos.update((n, (2, i)) for i, n in enumerate(Y))  # put nodes from Y at x=2
plt.figure(1,figsize=(50,50))
nx.draw(B, pos=pos,node_size=20)
plt.show()
