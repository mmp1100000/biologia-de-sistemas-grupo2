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

    """

        Esta función la vamos a utilizar para generar los diferentes gráficos de las redes. Para ello necesitamos un
        archivo con una estructura concreta. Este debe tener:

        Una columna para los farmacos y otra para los targets, ademas otra indicando si el target corresponde a cancer
        de mama. Que utilizaremos para plotear de diferentes colores. Y el p-valor que nos servirá para filtrar las
        interacciones.

    :param file: ruta del fichero

    :param type: selecciona el layout que prefieras:

        - spring
        - spectral
        - reingold
        - circular
        - kawai
        - random

    :param imageDim: Si tenemos muchas interacciones es recomendable aumentar la dimensión de la imagen. (Pasar un número a
    que hara que la imagen pase a ser del tamaño aXa)

    :return: networkx graph
    """

    # Leemos con pandas el fichero csv
    data = pd.read_csv(file, sep="\t")

    # Vemos el número de fármacos que tenemos
    unique_drugs = set(data['0'])
    total_number_of_tests = len(unique_drugs)
    df_node = pd.read_csv(file, sep="\t")

    # Filtramos el archivo por el p-valor
    data['p_value'] = pd.to_numeric(data['p_value'])
    df_node['p_value'] = pd.to_numeric(df_node['p_value'])
    data = data[data['p_value'] < (1 - conf_level)/total_number_of_tests]
    df_node = df_node[df_node['p_value'] < (1 - conf_level)/total_number_of_tests]

    warnings.filterwarnings('ignore')

    # Generamos un objeto gráfico donde añadiremos los nodos e interacciones
    G = nx.Graph(day="TodayNotTomorrow")

    # Contamos cuanto aparecen los fármacos y los targets para ajustar después el tamaño de los nodos
    veces = []
    veces1 = []
    for l in set(data['0'].tolist()):
        veces.append(data['0'].tolist().count(l))

    for l in set(data['1'].tolist()):
        veces1.append(data['1'].tolist().count(l) + 10)

    # Añadimos las interacciones
    df_edges = pd.DataFrame({
        'source': df_node['0'].tolist(),
        'target': df_node['1'].tolist(),
        'value': np.repeat(1, len(data['1'].tolist()))
    })

    #Aqui vamos a poner diferentes los colores
    vectorDif = data.sort_values('2', ascending=False).drop_duplicates('1', keep="first")
    colors_differ = vectorDif['2']
    groups = np.repeat(2, len(set(data['0'].values)),axis=0).tolist() + colors_differ.tolist()

    elements = list(set(data['0'])) + list(vectorDif['1'])

    # Añadimos los diferentes nodos
    df_nodes = pd.DataFrame({'name': elements,
                             'group': groups,
                             'nodesize': veces + veces1
                             })

    # Añadimos los nodos y conexiones
    for index, row in df_nodes.iterrows():
        G.add_node(row['name'], group=row['group'], nodesize=row['nodesize'])

    for index, row in df_edges.iterrows():
        G.add_weighted_edges_from([(row['source'], row['target'], row['value'])])

    # Seleccionamos los colores
    color_map = {0:'#DF0101' , 1: '#298A08', 2: '#0080FF'}

    # Ajustamos algunos parámetros para el layout
    plt.figure(figsize=(imageDim, imageDim))
    options = {
        'edge_color': '#FFDEA2',
        'width': 1,
        'with_labels': True,
        'font_weight': 'regular',
    }

    colors = [color_map[G.node[node]['group']] for node in G]
    sizes = [G.node[node]['nodesize'] * 40 for node in G]


    # Dependiendo del layout que queramos tendremos unos parámetros o otros
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
        pos_nodes = nx.spring_layout(G)
        nx.draw(G, node_color=colors, node_size=sizes, pos=pos_nodes, **options)
        node_attrs = nx.get_node_attributes(G, 'name')
        pos_attrs = {}
        
        for node, coords in pos_nodes.items():
            pos_attrs[node] = (coords[0], coords[1] + 0.08)
        
        nx.draw_networkx_labels(G, pos_attrs, labels=node_attrs)
        plt.savefig("outputs/net-kawai.png")
    else:
        nx.draw(G, node_color=colors, node_size=sizes, pos=nx.random_layout(G), **options)
        ax = plt.gca()
        ax.collections[0].set_edgecolor("#C0C0C0")
        plt.savefig("../outputs/net-random.png")


nnet_graph(inputfile, "kawai", imageDim=10)
# Si van a aparecer muchas conexiones, es mejor aumentar el tamaño de la ventana
# nnet_graph("../outputs/pvalued_interactions_0.95.tsv", "kawai", imageDim = 50)
