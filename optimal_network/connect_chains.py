import numpy as np
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from tqdm import tqdm
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import Voronoi
from itertools import combinations
from .utilities import toc, ILPNoSolutionError


def add_chains_to_strc(
        points: np.array,
        chains: np.array,
        chosen_trips: list,
        max_init_iter: int,
        debug: bool
    ):
    """
    Connects chains to a given structure by constructing optimal paths between fork nodes.
    The trips represent the structure as a tuple of v1, v2, c where v1, v2 are vertices and c is the chain

    Parameters
    ----------
    points : np.array
        Array of shape (n_points, 2) representing the coordinates of all points.
    chains : np.array
        Array assigning each point to a chain (same length as points).
    chosen_trips : list
        List of tuples (s, t, chain) specifying fork nodes (s, t) and their associated chain index.
    max_init_iter : int
        Maximum number of initialization iterations for chain optimization.
    debug : bool
        If True, show progress bars and debug information.

    Returns
    -------
    subgraph : networkx.Graph
        Graph with all chains optimally connected and node positions set as 'pos' attributes.
    """
    point_chain_df = _create_point_chain_df(points, chains, chosen_trips)
    subgraph = nx.Graph()
    progress = lambda x, **args: tqdm(x, **args) if debug else x
    for trip in progress(chosen_trips, desc="connect chains"):
        chain_edges = _connect_chain(point_chain_df, trip, max_init_iter)
        subgraph.add_edges_from(chain_edges)
    # add pos to subgraph
    for (node, row) in point_chain_df.iterrows():
        subgraph.nodes[node]["pos"] = row[["x", "y"]].values
    return subgraph


def _connect_chain(point_chain_df: pd.DataFrame, trip: tuple, max_init_iter: int):
    """Prepare data for a single chain and compute the optimal chain subgraph.

    This function extracts the points that belong to the requested chain (including
    the fork endpoints), tries to build a Delaunay-based sparse graph and solve
    the ILP on that graph. If the Delaunay graph is unsuitable or the ILP has no
    solution, it falls back to the complete graph for that chain.

    Parameters
    ----------
    point_chain_df : pd.DataFrame
        DataFrame with columns ['x', 'y', 'c', 'node', 'is_fork'] for all points.
    trip : tuple
        (s, t, chain) where s,t are fork node ids and chain is the chain index.
    max_init_iter : int
        Maximum initialization iterations for the ILP solver.

    Returns
    -------
    networkx.Graph.edges_view
        An edge-view (list-like) describing the edges selected for this chain.
    """
    # get the data of the chain
    s, t, chain = trip
    ends = [s, t]
    chain_data = point_chain_df.query('(c == @chain and (not is_fork)) or (node in @ends)')
    # try to construct the chain using delaunay_graph
    del_graph = _build_delaunay_graph(chain_data, ends)
    use_full_graph = del_graph is None
    if not use_full_graph:
        try:
            chain = _optimal_chain(del_graph, s, t, max_init_iter=max_init_iter)
        except ILPNoSolutionError as e_dt:
            use_full_graph = True
    if use_full_graph:
        # fallback: use the complete graph
        full_graph = _build_full_graph(chain_data, ends)
        chain = _optimal_chain(full_graph, s, t, max_init_iter=max_init_iter)
    return chain.edges


def key(e):
    """Normalized key for an edge tuple: (min, max).

    Useful for canonical comparisons of undirected edges.
    """
    return tuple(sorted(e))

def key_edges(edges):
    """Return a canonical tuple representation for a collection of edges.

    Each edge is normalized with :func:`key` and the list is sorted so it can be
    used as a hashable key (for example when tracking added constraints).
    """
    return tuple(sorted([key(e) for e in edges]))



def _create_point_chain_df(points: np.array, chains: np.array, chosen_trips: list) -> pd.DataFrame:
    """Create a DataFrame describing points, their chain assignment and fork flags.

    The returned DataFrame contains columns ['x','y','c','node','is_fork'] where
    'is_fork' marks nodes that are endpoints (forks) referenced in chosen_trips.
    """
    point_chain_df = pd.DataFrame(np.c_[points, chains], columns=["x", "y", "c"])
    point_chain_df['node'] = list(range(len(points)))
    point_chain_df['node'] = point_chain_df['node'].apply(int).astype('int')
    all_forks = set([fork for trip in chosen_trips for fork in trip[:2]])
    point_chain_df['is_fork'] = point_chain_df['node'].isin(all_forks)
    return point_chain_df


def _build_delaunay_graph(point_chain_df: pd.DataFrame, ends: list) -> nx.Graph:
    """Build a sparse geometric graph using Voronoi/Delaunay ridge connections.

    Returns a weighted NetworkX graph connecting nearby chain points; returns
    ``None`` if Voronoi/Delaunay computation fails.
    """
    points = point_chain_df[["x", "y"]].values
    nodes = [int(x) for x in point_chain_df["node"]]
    dist = squareform(pdist(points))
    try:
        vor = Voronoi(points)
    except:
        return None
    ends_idx = [nodes.index(end) for end in ends]
    edges = [(nodes[u], nodes[v], dist[u, v]) for u, v in vor.ridge_points]
    edges += [(nodes[u], nodes[end], dist[u, end]) for end in ends_idx for u in range(len(nodes))]
    del_graph = nx.Graph()
    del_graph.add_weighted_edges_from(edges)
    return del_graph

def _build_full_graph(data: pd.DataFrame, ends: list) -> nx.Graph:
    """Build a complete weighted graph between the provided points.

    This is used as a fallback when the Delaunay-based sparse graph is not available.
    """
    points = data[["x", "y"]].values
    nodes = [int(x) for x in data["node"]]
    dist = squareform(pdist(points))
    edges = [(nodes[u], nodes[v], dist[u, v]) for u, v in combinations(range(len(nodes)), 2)]
    full_graph = nx.Graph()
    full_graph.add_weighted_edges_from(edges)
    return full_graph

def _optimal_chain(G: nx.Graph, s, t, max_init_iter: int=2):
    """High-level wrapper that computes an optimal s--t chain in graph G.

    Attempts to solve using iterative ILP calls while preferring previously found
    edges to improve connectivity in repeated solves.
    """
    prefer_prev_edges = max_init_iter is not None
    subgraph = _solve_optimal_chain(
        G, s, t,
        max_init_iter=max_init_iter,
        prefer_prev_edges=prefer_prev_edges,
        start_edges=None
    )
    if prefer_prev_edges:
        comps = list(nx.k_edge_components(subgraph, 2))
        if len(comps) == 1:
            return subgraph
        subgraph = _solve_optimal_chain(
            G, s, t,
            max_init_iter=None,
            prefer_prev_edges=True,
            start_edges=subgraph.edges
        )
    return subgraph


def _solve_optimal_chain(
        G: nx.Graph, s, t,
        max_init_iter: int=2,
        prefer_prev_edges: bool=False,
        start_edges=None
    ):
    """Solve ILP to find a minimum-weight s--t subgraph satisfying degree constraints.

    Implements the ILP model using Gurobi binary variables for edges and iteratively
    adds connectivity cuts until the selected edges form a connected s--t chain.

    Parameters
    ----------
    G : nx.Graph
        Weighted input graph (edge attribute 'weight' required).
    s, t : hashable
        Source and sink node ids within G.
    max_init_iter : int or None
        Maximum number of outer iterations for adding connectivity cuts. If None,
        performs a single solve with stronger preference for previous edges.
    prefer_prev_edges : bool
        When True, previously preferred edges are excluded from the objective to
        bias the solver towards keeping them.
    start_edges : iterable or None
        Optional list of edges to warm-start or prefer.

    Returns
    -------
    networkx.Graph
        Edge-subgraph of G containing the selected edges.
    """
    # solve for a graph with len 0
    res = _optimal_chain_trivial(G, s, t)
    if res is not None:
        return res
    
    # init model
    m = gp.Model("connect_chain")
    m.Params.OutputFlag = 0
    constrain = set()

    if start_edges:
        prefer_edges = list(map(key, start_edges))
    else:
        prefer_edges = []

    edges = list(map(key, G.edges))
    x = m.addVars(edges, vtype=GRB.BINARY, name="x")  # edges

    # n_nodes - 1 edges
    n_nodes = len(G.nodes) - 1 if s != t else len(G.nodes)
    m.addConstr(gp.quicksum(x[e] for e in edges) == n_nodes, name="edge_count")

    # degree constrained
    for v in G.nodes:
        if s == t:
            degree = 2
        else:
            degree = 1 if v in [s, t] else 2
        m.addConstr(gp.quicksum(x[key(e)] for e in G.edges(v)) == degree, name="degree")

    # objective
    m.setObjective(gp.quicksum(G.edges[e]['weight'] * x[e] for e in edges if e not in prefer_edges), GRB.MINIMIZE)
    m.optimize()

    # run model
    if max_init_iter is None:
        max_init_iter = 1e7
    it = 0
    while it < max_init_iter:
        status = m.Status
        if status != GRB.Status.OPTIMAL:
            raise ILPNoSolutionError(f"no solution found for chain for chain {s, t}")
        
        # add connectivity constrain
        selected_edges = [e for e in edges if x[e].X > 0.5]
        subgraph = nx.from_edgelist(selected_edges)
        if prefer_prev_edges:
            prefer_edges = selected_edges
        comps = list(nx.connected_components(subgraph))
        if len(comps) == 1:
            break
        for comp in comps:
            cut_edges = [e for e  in edges if len(set(e) & comp) == 1]
            if key_edges(cut_edges) not in constrain:
                m.addConstr(gp.quicksum(x[e] for e in cut_edges) >= 1)
                constrain.add(key_edges(cut_edges))
        m.setObjective(gp.quicksum(G.edges[e]['weight'] * x[e] for e in edges if e not in prefer_edges), GRB.MINIMIZE)
        m.optimize()
    
    # return subgraph
    selected_edges = [e for e in edges if x[e].X > 0.5]
    subgraph = nx.edge_subgraph(G, selected_edges).copy()
    return subgraph



def _optimal_chain_trivial(G: nx.Graph, s, t):
    """Quick trivial-case check for tiny graphs.

    Returns an immediate solution when G is empty (connect s--t directly) or
    ``None`` when no trivial solution exists.
    """
    if len(G) == 0:
        if s != t:
            return nx.from_edgelist([(s, t)])
        return nx.Graph()
    return None