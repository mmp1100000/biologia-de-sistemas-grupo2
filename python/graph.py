import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

data = pd.read_csv('../outputs/interactions.tsv', sep="\t", header=None)

import warnings
warnings.filterwarnings('ignore')

G = nx.Graph(day="Stackoverflow")

veces = []
veces1 = []
for l in set(data[0].tolist()):
    veces.append(data[0].tolist().count(l))

for l in set(data[1].tolist()):
    veces1.append(data[1].tolist().count(l))


df_node = pd.read_csv('../outputs/interactions.tsv', sep="\t", header=None)

df_edges = pd.DataFrame({
    'source': df_node[0].tolist(),
    'target': df_node[1].tolist(),
    'value': np.repeat(1, len(data[1].tolist()))
})
elements = list(set(data[0])) + list(set(data[1]))
groups = np.repeat(0, len(set(data[0].values)), axis=0).tolist() + np.repeat(1, len(set(data[1].values)), axis=0).tolist()
df_nodes = pd.DataFrame({'name': elements ,
                        'group':  groups,
                         'nodesize': veces + veces1
                         })
for index, row in df_nodes.iterrows():
    G.add_node(row['name'], group=row['group'], nodesize=row['nodesize'])


for index, row in df_edges.iterrows():
    G.add_weighted_edges_from([(row['source'], row['target'], row['value'])])

color_map = {0: '#f09494', 1: '#d2f5f0'}

plt.figure(figsize=(50, 50))
options = {
    'edge_color': '#FFDEA2',
    'width': 1,
    'with_labels': True,
    'font_weight': 'regular',
}


colors = [color_map[G.node[node]['group']] for node in G]
sizes = [G.node[node]['nodesize']*10 for node in G]


"""

Using the spring layout : 
- k controls the distance between the nodes and varies between 0 and 1
- iterations is the number of times simulated annealing is run
default k=0.1 and iterations=50

"""

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.spring_layout(G, k=0.25, iterations=50), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-spring.png")

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.spectral_layout(G), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-spectral.png")

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.fruchterman_reingold_layout(G), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-reingold.png")

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.random_layout(G), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-random.png")

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.circular_layout(G), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-circular.png")

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.shell_layout(G), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-shell.png")

nx.draw(G, node_color=colors, node_size=sizes, pos=nx.kamada_kawai_layout(G), **options)
ax = plt.gca()
ax.collections[0].set_edgecolor("#C0C0C0")
plt.savefig("../outputs/NNETGraph/net-kawai.png")