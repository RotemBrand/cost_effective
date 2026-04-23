import numpy as np
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from .cluster_chains import balanced_kmeans_gurobi
from .construct_strc import get_optimal_strc_trips
from .connect_chains import add_chains_to_strc
from .utilities import toc

def optimal_network_from_points(
    points: np.array,
    r: int,
    source_node: int=None,
    chain_len_sigma: float=0,
    kmeans_max_iter: int=7,
    strc_n_init_iters: int=2,
    strc_exact_vertices: bool=False,
    chain_n_init_iters: int=2,
    seed: int=7,
    debug: bool=False
):
    r"""
    Constructs an optimal network graph from a set of points with low total length.
    An `optimal graph` is 2 connected with almost equal chains and a structure that it
    3 connected

    Parameters
    ----------
    points : np.array
        Array of shape (n_points, 2) representing the coordinates of the points.
    r : int
        Redundancy parameter; controls the number of chains and structure nodes.
    chain_len_sigma : float, optional
        Standard deviation for chain length regularization in clustering (default: 0).
    kmeans_max_iter : int, optional
        Maximum number of iterations for k-means clustering (default: 7).
    strc_n_init_iters : int, optional
        Number of initialization iterations for structure optimization -
        set a low value to increase preformance over weight (default: 2).
    strc_exact_vertices : bool, optional
        If True, use exact vertices for structure optimization, otherwise, use all the vertices of the delauaney triangulation (default: False).
    chain_n_init_iters : int, optional
        Number of initialization iterations for chain construction-
        set a low value to increase preformance over weight (default: 2).
    seed : int, optional
        Random seed for reproducibility (default: 7).
    debug : bool, optional
        If True, print timing and debug information (default: False).

    Returns
    -------
    subgraph : networkx.Graph
        The constructed optimal network as a NetworkX graph.
    """
    # init
    if r > 1:
        n_nodes = 2 * (r - 1)
        n_edges = 3 * (r - 1)
    elif r == 1:
        n_nodes, n_edges = 1, 1
    else:
        raise ValueError("r should be more than 0")
    # clusters the points into chains
    toc(debug=debug)
    chains, centers = balanced_kmeans_gurobi(
        X=points,
        k=n_edges,
        max_iter=kmeans_max_iter,
        chain_len_sigma=chain_len_sigma,
        random_state=seed,
    )
    toc("balanced_kmeans_gurobi", debug=debug)
    # build structure
    chosen_trips = get_optimal_strc_trips(
        centers=centers,
        points=points,
        n_nodes=n_nodes,
        strc_n_init_iters=strc_n_init_iters,
        exact_vertices=strc_exact_vertices,
        source_node=source_node,
        debug=debug
    )
    toc("get_optimal_strc_trips", debug=debug)
    # add chains
    subgraph = add_chains_to_strc(
        points=points,
        chains=chains,
        chosen_trips=chosen_trips,
        max_init_iter=chain_n_init_iters,
        debug=debug
    )
    toc("add_chains_to_strc", debug=debug)
    return subgraph
