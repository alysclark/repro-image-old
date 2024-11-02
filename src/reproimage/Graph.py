import copy
import networkx as nx
import numpy as np
import SimpleITK as sitk
from utils import cartesian_product
from numba import jit
from scipy import spatial
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

def get_junction_nodes(graph):
    """
    :param graph: nx.graph
    :return: a tuple containing the number of junction nodes, and the set of junction nodes
    Junction nodes are any node with a degree higher than 2
    """
    deg = nx.degree(graph)
    degree_node_set = []
    for node in graph.nodes:
        if deg[node] > 2:
            degree_node_set.append(node)
    return len(degree_node_set), degree_node_set

def largest_ccmp_nx(graph):
    """
    :param graph:
    :return: The largest connected component of the input graph
    """
    largest_cc = max(nx.connected_components(graph), key=len)
    S = graph.subgraph(largest_cc).copy()
    return S

def give_me_tree(J_graph, full_graph, coords, segmentation_image):
    """This has been tried with just eccentricity but that worked like  so now we will try to incorporate
     length, and if that doesn't work then I will use average diameter, this function assumes the graph has a mapped
     junction node density field"""
    junction_density_field, junction_nodes = calculate_junction_node_density_field(full_graph, coords, radius=35)
    map_node_field_graph(J_graph, junction_density_field, junction_nodes, label='JDF')
    norm_vector_perms = np.array([0, -0.5, 0.5])
    perm_array_for_normals = cartesian_product(norm_vector_perms, norm_vector_perms)
    perm_array_for_normals = np.delete(perm_array_for_normals, 0, 0)
    seg_array = sitk.GetArrayFromImage(segmentation_image)
    cycle_basis = nx.cycle_basis(J_graph)
    removed_edges = []
    for cycle in cycle_basis:
        edge_eccen = []
        edge_len = []
        edge_rad = []
        edge_JDF = []
        edges = []
        edge_chorionic_maternal = []
        cycle_length = len(cycle)
        # iterate over edges in cycle
        for ind in range(0, len(cycle)):
            u, v = (cycle[ind], cycle[(ind + 1) % cycle_length])
            if J_graph.get_edge_data(u, v) == None or J_graph.get_edge_data(u, v) == {}:
                edge_nodes = [u, v]
            else:
                # print(J_graph.get_edge_data(u, v))
                edge_nodes = J_graph.get_edge_data(u, v)['path']
            if len(edge_nodes) == 1:
                edge_nodes = [u, edge_nodes[0], v]
            # print(f"u: {u}: {J_graph.degree[u]}, v: {v}: {J_graph.degree[v]}")
            edge_node_coords = coords[edge_nodes]
            edge_chorionic_maternal.append(edge_node_coords.max(axis=0)[1])
            edge_JDF.append(max([J_graph.nodes[x]['JDF'] for x in [u,v]]))
            edge_len.append(len(edge_nodes))
            edge_eccen.append(eccentricty_path(edge_node_coords, perm_array_for_normals, seg_array))
            edge_rad.append(radius_path(edge_node_coords, perm_array_for_normals, seg_array))
            edges.append((u, v))
        normalised_rad = [1 - x/max(edge_rad) for x in edge_rad]
        normalised_length = [1-x/max(edge_len) for x in edge_len]
        normalised_JDF = [x/max(edge_JDF) for x in edge_len]
        normalised_maternal_chorionic = [x/max(edge_chorionic_maternal) for x in edge_chorionic_maternal]
        eccen_length_product = []
        for norm_JDF, eccen, norm_rad, norm_len, norm_matern in zip(normalised_JDF, edge_eccen, normalised_rad, normalised_length, normalised_maternal_chorionic):
            eccen_length_product.append(eccen*norm_matern)

        max_eccen_index = eccen_length_product.index(max(eccen_length_product))
        u_remove, v_remove = edges[max_eccen_index]
        removed_edges.append((u_remove, v_remove))
    J_graph.remove_edges_from(removed_edges)
    return J_graph, removed_edges

def eccentricty_path(pixel_coords, perm_array, seg_img):
    """
    :param pixel_coords:
    :param perm_array:
    :param seg_img:
    :return: the mean eccentricity across the path defined by the pixel coordinate inputs
    """
    eccen = []
    edge_measure = np.zeros((pixel_coords.shape[0]-1, 8))
    for p,_ in enumerate(pixel_coords[:-1]):
        edge_measure = find_distances_using_normal(pixel_coords[p], pixel_coords[p+1],seg_img,
                                               perm_array)
        e_max = edge_measure.max(axis=0)
        e_min = edge_measure.min(axis=0)
        if np.isnan(e_min):
            e_min = 1.0
        eccen.append(e_max/e_min)
    return np.mean(eccen)

def radius_path(pixel_coords, perm_array, seg_img, metric='mean'):
    """
    :param pixel_coords:
    :param perm_array:
    :param seg_img:
    :param metric:
    :return: returns a measure of the path radius determined by given metric
    """

    radius = []
    for p,_ in enumerate(pixel_coords[:-1]):
        if np.array([x==z for x,z in zip(pixel_coords[p], pixel_coords[p+1])]).all():
            edge_measure = 0
        else:
            edge_measure = find_distances_using_normal(pixel_coords[p], pixel_coords[p+1],seg_img, perm_array)
        radius.append(np.mean(edge_measure))
    if metric == 'mean':
        return np.mean(radius)
    elif metric == 'max':
        return max(radius)

@jit(nopython=True)
def find_distances_using_normal(coord1, coord2, VolumeImage, perm_array):

    # get centre line vector
    centre = (coord1 - coord2).astype(np.double)/ np.linalg.norm((coord1 - coord2).astype(np.double))
    numSamples = perm_array.shape[0]
    distances = np.ones(numSamples)
    normal = np.zeros(3)
    for i in range(0,numSamples):
        # Randomly assign normal vector, using the dot product rule (centre.normal==0) and avoiding div0 errors
        if centre[2]!= 0:
            normal[0] = perm_array[i,0]
            normal[1] = perm_array[i,1]
            normal[2] = -(centre[0] * normal[0] + centre[1] * normal[1])/ centre[2]

        elif centre[0] != 0:
            normal[2] = perm_array[i,0]
            normal[1] = perm_array[i,1]
            normal[0] = -(centre[2] * normal[2] + centre[1] * normal[1]) / centre[0]

        else: # centre[1]!= 0:
            normal[0] = perm_array[i,0]
            normal[2] = perm_array[i,1]
            normal[1] = -(centre[0] * normal[0] + centre[2] * normal[2]) / centre[1]

        normal = normal / np.linalg.norm(normal)

        # Find distances
        step = 0
        counter = 0
        currentValue = 1
        while (currentValue == 1) & (counter < 200): # check if in vessel (plus arbitrary check)

             step = step + 0.1 # step update by 1/5 of a voxel (could increase in order to speed up)
             counter = counter + 1
             y = (coord1).astype(np.double) + step*normal # take step in direction of normal vector
             currentPosition = np.empty_like(y)
             np.round_(y, 0, currentPosition)
             if(int(currentPosition[0])>np.shape(VolumeImage)[0]-1):
                 currentValue = 0
             elif(int(currentPosition[1])>np.shape(VolumeImage)[1]-1):
                 currentValue = 0
             elif(int(currentPosition[2])>np.shape(VolumeImage)[2]-1):
                 currentValue = 0
             else:
                 currentValue = VolumeImage[int(currentPosition[0]), int(currentPosition[1]), int(currentPosition[2])]
        distances[i] = step - 0.1
    distances[distances == 0] = 0.5 # set the radius for points with a zero radius to half a voxel, as that is the
    # theoretical minimum radius
    return distances

def map_node_field_graph(graph, field, nodes, label='field'):
    for node, val in zip(nodes, field):
        if graph.has_node(node):
            nx.set_node_attributes(graph, {node: val}, name=label)

def calculate_junction_node_density_field(graph, coords, radius=35):
    n_junction, junction_nodes = get_junction_nodes(graph)
    junc_coords = coords[junction_nodes]
    tree = spatial.KDTree(junc_coords)
    results = tree.query_ball_point([x for x in junc_coords], radius)
    density_values = [len(x) for x in results]
    return density_values, junction_nodes