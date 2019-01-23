import sys

import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import warnings

from plotly.offline.offline import matplotlib

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


def bipartite_graph(file):
    # plot dimensions
    matplotlib.rcParams['figure.figsize'] = [15, 20]  # for square canvas
    matplotlib.rcParams['figure.subplot.left'] = 0
    matplotlib.rcParams['figure.subplot.bottom'] = 0
    matplotlib.rcParams['figure.subplot.right'] = 1
    matplotlib.rcParams['figure.subplot.top'] = 1

    # read flagged interactions: colums 0 and 1 as nodes, rows as interactions
    data = pd.read_csv(file, sep="\t")
    interactions = [tuple([row[1], row[2], row[4]]) for row in data.values.tolist() if row[3] == 1]
    drug_pval = set([tuple([row[0], row[2]]) for row in interactions])
    drug_pval = [row[1] for row in drug_pval]
    drug_nodes = list(set([row[0] for row in interactions]))
    gene_nodes = list(set([row[1] for row in interactions]))
    interactions = [tuple([row[0], row[1]]) for row in interactions]
    print(len(interactions))
    print(len(drug_nodes))
    print(len(gene_nodes))

    # create bipartite graph
    G = nx.Graph()
    # Add nodes with the node attribute "bipartite"
    G.add_nodes_from(drug_nodes, bipartite=0)
    G.add_nodes_from(gene_nodes, bipartite=1)
    # Add edges only between nodes of opposite node sets
    G.add_edges_from(interactions)
    nx.is_connected(G)

    # distribute node positions according to their bipartite 'class' 0-drugs, 1-genes
    X, Y = nx.bipartite.sets(G)
    pos = dict()
    posX = dict()
    posY = dict()
    posX.update((n, (0, i * 10)) for i, n in enumerate(X))
    posY.update((n, (0.5, i * 100)) for i, n in enumerate(Y))
    pos.update((n, (0, i * 10)) for i, n in enumerate(X))
    pos.update((n, (0.5, i * 100)) for i, n in enumerate(Y))

    bipartite_colors = np.array(list(reversed(sorted(list(nx.bipartite.color(G).values()), key=int))))
    print(bipartite_colors)

    # draw in the positions
    nx.draw(G, pos=pos, with_labels=True,
            edge_color="Green",
            edge_cmap=plt.get_cmap('Greens'),
            node_color="Red",
            cmap=plt.get_cmap('Reds'))

    # nx.draw_networkx_nodes(X, pos=posX, with_labels=True, node_color=np.array(drug_pval), cmap=plt.get_cmap('Greens'))
    # nx.draw_networkx_nodes(Y, pos=posY, with_labels=True, node_color="Blue")
    # nx.draw_networkx_labels(G, pos=pos)
    # nx.draw_networkx_edges(G, pos=pos)

    plt.title('All breast-cancer-genes interactions with GSK drugs')
    plt.savefig("../outputs/NNETGraph/bipartite.png")


nnet_graph(inputfile, "kawai", imageDim=10)
bipartite_graph(inputfile)

# Si van a aparecer muchas conexiones, es mejor aumentar el tamaño de la ventana
# nnet_graph("../outputs/pvalued_interactions_0.95.tsv", "kawai", imageDim = 50)
