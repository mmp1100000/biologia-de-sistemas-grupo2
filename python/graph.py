import sys

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings

if len(sys.argv) < 3:
    exit("pvalued_interactions tsv filepath and confidence interval index required.")
inputfile = sys.argv[1]
conf_level = float(sys.argv[2])


def nnet_graph(file, type, imageDim=20):
    data = pd.read_csv(file, sep="\t")
    df_node = pd.read_csv(file, sep="\t")
    data['p_value'] = pd.to_numeric(data['p_value'])
    df_node['p_value'] = pd.to_numeric(df_node['p_value'])

    data = data[data['p_value'] < (1 - conf_level)]
    df_node = df_node[df_node['p_value'] < (1 - conf_level)]

    warnings.filterwarnings('ignore')

    G = nx.Graph(day="TodayNotTomorrow")

    veces = []
    veces1 = []
    for l in set(data['0'].tolist()):
        veces.append(data['0'].tolist().count(l))

    for l in set(data['1'].tolist()):
        veces1.append(data['1'].tolist().count(l))

    df_edges = pd.DataFrame({
        'source': df_node['0'].tolist(),
        'target': df_node['1'].tolist(),
        'value': np.repeat(1, len(data['1'].tolist()))
    })
    elements = list(set(data['0'])) + list(set(data['1']))
    groups = np.repeat(0, len(set(data['0'].values)), axis=0).tolist() + np.repeat(1, len(set(data['1'].values)),
                                                                                   axis=0).tolist()
    df_nodes = pd.DataFrame({'name': elements,
                             'group': groups,
                             'nodesize': veces + veces1
                             })
    for index, row in df_nodes.iterrows():
        G.add_node(row['name'], group=row['group'], nodesize=row['nodesize'])

    for index, row in df_edges.iterrows():
        G.add_weighted_edges_from([(row['source'], row['target'], row['value'])])

    color_map = {0: '#f09494', 1: '#d2f5f0'}

    plt.figure(figsize=(imageDim, imageDim))
    options = {
        'edge_color': '#FFDEA2',
        'width': 1,
        'with_labels': True,
        'font_weight': 'regular',
    }

    colors = [color_map[G.node[node]['group']] for node in G]
    sizes = [G.node[node]['nodesize'] * 40 for node in G]

    if type == "spring":
        """

        Using the spring layout : 
        - k controls the distance between the nodes and varies between 0 and 1
        - iterations is the number of times simulated annealing is run
        default k=0.1 and iterations=50

        """
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.spring_layout(G, k=0.25, iterations=50), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-spring.png")

    elif type == "spectral":
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.spectral_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-spectral.png")
    elif type == "reingold":
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.fruchterman_reingold_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-reingold.png")
    elif type == "circular":
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.circular_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-circular.png")
    elif type == "shell":
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.shell_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-shell.png")
    elif type == "kawai":
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.kamada_kawai_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("outputs/net-kawai.png")
    else:
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.random_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-random.png")


nnet_graph(inputfile, "kawai", imageDim=10)
# Si van a aparecer muchas conexiones, es mejor aumentar el tamaño de la ventana
# nnet_graph("../outputs/pvalued_interactions_0.95.tsv", "kawai", imageDim = 50)
