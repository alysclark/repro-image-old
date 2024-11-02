import numpy as np

def remap_node_field_for_vis(graph, field):
    """
    graph - nx.Graph
    field - np.array
    assuming field corresponds to graph.nodes ordering
    """
    graph_node_mapping = {k:v for k, v in zip(graph.nodes, range(graph.number_of_nodes()))}
    nodes = np.unique(np.array(graph.edges()).flatten())
    new_field = []
    for node in nodes:
        new_field.append(field[graph_node_mapping[node]])
    return new_field


def generate_visualisation_arrays(coords, edges):
    """coords - coordinate array for glaboal node numbers
    edges, node connections for global node numbers"""
    nodes = np.unique(edges.flatten())
    nodes.sort()
    new_coords = coords[nodes]
    new_map = {p: c for c, p in enumerate(nodes)}
    remapped_connections = []
    for element in edges:
        temp = []
        for node in element:
            temp.append(new_map[node])
        remapped_connections.append(temp)
    remapped_connections = np.array(remapped_connections)
    padding = np.empty(remapped_connections.shape[0], int) * 2
    padding[:] = 2
    connections_with_padding = np.vstack((padding, remapped_connections.T)).T
    return new_coords, connections_with_padding