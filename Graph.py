import copy
import networkx as nx

def get_node_degree_set(graph, degreeValue):
    """
    :param graph: undirected graph
    :param degreeValue:integer defining the degree of each node in the output set
    :return: a set of nodes in the input graph that have the degree value specified
    """
    deg = nx.degree(graph)
    degree_node_set = []
    for node in graph.nodes:
        if deg[node] == degreeValue:
            degree_node_set.append(node)
    return degree_node_set

def junction_node_subgraph(graph):
    """
    :param graph:
    :return: a series reduced tree with no degree 2 nodes
    """
    jgraph = copy.deepcopy(graph)
    # Iterate through the edges and remove attributes
    for u, v, attrs in jgraph.edges(data=True):
        for attr_key in list(attrs.keys()):
            jgraph[u][v].pop(attr_key)

    parallel_edges = True
    while parallel_edges:
        parallel_edges = False
        for node in graph.nodes():
            if jgraph.has_node(node):
                if graph.degree[node] == 2: # for degree 2 nodes

                    new_edge = []
                    # assert len(list(jgraph.neighbors((item[0])))) == 2
                    for n_node in jgraph.neighbors(node):
                        new_edge.append(n_node)
                    # print(f"Processing node: {node}, with neighbors: {new_edge}")
                    if len(new_edge) == 2:
                        new_path = []
                        for edge in jgraph.edges(node):
                            if jgraph.get_edge_data(edge[0], edge[1]) != {}:
                                # print(jgraph.get_edge_data(edge[0], edge[1]))
                                new_path = jgraph.get_edge_data(edge[0], edge[1])['path']
                        jgraph.remove_node(node)
                        u, v = new_edge
                        # print(f'adding edge: {u,v}')
                        if not jgraph.has_edge(u, v):
                            jgraph.add_edge(u, v)
                            if jgraph.get_edge_data(u, v) == {}:
                                new_path.append(node)
                                nx.set_edge_attributes(jgraph, {(u, v): new_path}, name='path')
                        else:
                            parallel_edges = True
                    elif len(new_edge) == 1:
                        for edge_node in [node, new_edge[0]]:
                            if graph.degree[edge_node] == 2:
                                jgraph.remove_node(edge_node)
                #else:
                    # print(f"not processing node: {node}, with degree: {jgraph.degree[node]}, and neighbors: {list(jgraph.neighbors(node))}")
    return jgraph