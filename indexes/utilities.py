from __future__ import annotations

import networkx as nx
import numpy as np
from itertools import combinations
import indexes.graph_rel as GR
from typing import Union, Dict, List, Optional, Tuple, Literal



###### skeleton functions ######
def get_skeleton_graph(graph, sources: Optional[List]=None):
    if sources is None:
        sources = [defult_sources(graph)]
    # validate sources
    for source in sources:
        if source not in graph.nodes:
            raise ValueError(f"Source {source} not found in graph")    
    # init skeleton
    skeleton = nx.MultiGraph()
    graph_copy = nx.MultiGraph(graph)
    unvisited_forks = set([node for node, d in graph.degree if d != 2] + sources)
    # make DFS starting from one of the forks to add a chain
    while unvisited_forks:
        start_node = next(iter(unvisited_forks))
        curr = start_node
        edges = []
        first_iter = True
        # dfs
        while first_iter or (curr not in unvisited_forks):
            if graph_copy.degree(curr) == 0:
                if curr in unvisited_forks:
                    unvisited_forks.remove(curr)
                break
            prev, curr, k = next(iter(graph_copy.edges(curr, keys=True)))
            edges.append((prev, curr, k))
            graph_copy.remove_edge(prev, curr, key=k)
            first_iter = False          
        # add the new chain to skeleton 
        if not first_iter:
            chain_attr = _get_chain_attributes(graph, edges)
            skeleton.add_edge(start_node, curr, **chain_attr)
    # add forks attributes
    for fork in skeleton.nodes:
        skeleton.nodes[fork].update(graph.nodes[fork])


    return skeleton

def defult_sources(graph):
    try:
        source = next(node for node, d in graph.degree if d >= 3)
    except StopIteration:
        source = next(iter(graph.nodes))
    return source

def _get_chain_attributes(graph, edges: list) -> dict:
    # edges
    if isinstance(graph, nx.MultiGraph):
        adjusted_edges = edges
    elif isinstance(graph, nx.Graph):
        adjusted_edges = [e[:2] for e in edges] 
    else:
        raise ValueError(f"graph type sould be nx.Graph or nx.MultiGraph not {type(graph)}")
    # nodes
    nodes = [e[0] for e in adjusted_edges]
    nodes.append(adjusted_edges[-1][1])
    # length
    length = len(adjusted_edges)
    # subgraph
    subgraph = nx.edge_subgraph(graph, adjusted_edges).copy()
    return {
        'edges': adjusted_edges,
        'nodes': nodes,
        'length': length,
        'subgraph': subgraph
    }

# old skeleton function

def create_gr_of_the_skeleton_graph(gr):
    skeleton = get_skeleton_graph(gr.graph, [gr.source])
    # calculate the failing prob of each edge 1 - prod(1 - p_e for e in edge.subgraph.edges)
    probs = {
        edge: -np.prod(
            -np.array(list(nx.get_edge_attributes(edge_data['subgraph'], 'prob').values())) + 1,
            axis=0
        ) + 1
        for edge, edge_data in skeleton.edges.items()
    }
    gr_skeleton = GR.GraphRel(
        graph=skeleton,
        nodes_weight={node: gr.graph.nodes[node]["weight"] for node in skeleton.nodes},
        edges_prob=probs,
        max_fail=gr.max_fail,
        sources=[gr.source],
    )
    return gr_skeleton



###### prob by length ######
def edge_probs_by_length(
    graph: nx.Graph,
    *,
    p: float,
    prob_attr: str | None = None,
    pos_attr: str = "pos",
    mode: Literal["rate", "mean"] = "rate",
    tol: float = 1e-5,
    max_iter: int = 100,
) -> Tuple[Dict[Tuple[object, object], float], float]:
    """
    Compute per-edge failure probabilities based on geometric edge length.

    For edge e with length l_e:
        P(fail) = 1 - exp(-l_e * rate)

    Parameters
    ----------
    graph : nx.Graph
    p : float
        If mode="rate": p is the failure rate per unit length (>= 0).
        If mode="mean": p is the target mean edge failure probability in (0, 1),
                        and we solve for 'rate' by bisection.
    prob_attr : str | None
        If provided, also writes the computed probability to graph[u][v][prob_attr].
    pos_attr : str
        Node attribute containing position as (x, y).
    mode : {"rate", "mean"}
    tol, max_iter : solver settings for mode="mean"

    Returns
    -------
    probs : dict
        { (min(u,v), max(u,v)) : p_e }  with sorted tuple edge keys.
    rate: float
        the failing rate per length used to calculate
    """
    if graph.number_of_edges() == 0:
        return {}

    pos = nx.get_node_attributes(graph, pos_attr)
    if len(pos) != graph.number_of_nodes():
        missing = [n for n in graph.nodes if n not in pos][:5]
        raise KeyError(f"Missing '{pos_attr}' for some nodes, e.g. {missing}")

    # Stable edge list + keys
    edges = list(graph.edges())
    keys = [tuple(sorted((u, v))) for (u, v) in edges]  # assumes nodes are orderable

    # Vectorized lengths
    xy_u = np.array([pos[u] for (u, _) in edges], dtype=float)
    xy_v = np.array([pos[v] for (_, v) in edges], dtype=float)
    lengths = np.hypot(xy_v[:, 0] - xy_u[:, 0], xy_v[:, 1] - xy_u[:, 1])

    if mode == "rate":
        rate = float(p)
        if rate < 0:
            raise ValueError("In mode='rate', p must be >= 0 (rate per unit length).")
    elif mode == "mean":
        rate = _solve_rate_for_mean_failure_prob(lengths, float(p), tol=tol, max_iter=max_iter)
    else:
        raise ValueError("mode must be 'rate' or 'mean'.")

    probs_arr = 1.0 - np.exp(-lengths * rate)
    probs_arr = np.clip(probs_arr, 0.0, 1.0)

    probs = {k: float(pe) for k, pe in zip(keys, probs_arr)}

    if prob_attr is not None:
        for (u, v), pe in zip(edges, probs_arr):
            graph[u][v][prob_attr] = float(pe)
    return probs, rate

def _solve_rate_for_mean_failure_prob(
    lengths: np.ndarray,
    target_p: float,
    *,
    tol: float = 1e-5,
    max_iter: int = 100,
) -> float:
    """
    Find rate >= 0 such that mean(1 - exp(-lengths * rate)) == target_p,
    using bisection (monotone in rate).

    lengths: 1D array of nonnegative edge lengths.
    target_p: desired mean failure probability in (0, 1).
    """
    lengths = np.asarray(lengths, dtype=float)
    if lengths.ndim != 1:
        raise ValueError("lengths must be a 1D array")
    if np.any(lengths < 0):
        raise ValueError("lengths must be nonnegative")
    if target_p == 0.0:
        return 0.0
    if target_p == 1.0:
        raise ValueError("target_p must be in [0, 1)")

    if lengths.size == 0 or float(lengths.max()) == 0.0:
        # no length (or all zero-length edges) => failure prob always 0
        return 0.0

    # mean_prob(rate) = mean(1 - exp(-l*rate))
    def mean_prob(rate: float) -> float:
        return float(np.mean(1.0 - np.exp(-lengths * rate)))

    # Small-rate approximation: 1-exp(-l r) ~ l r
    # mean ~ r * mean(l)  => r ~ target_p / mean(l)
    mean_l = float(np.mean(lengths))
    guess = target_p / mean_l

    lo = 0
    hi = max(guess, 1e-16)

    # Ensure bracket: mean_prob(lo) <= target <= mean_prob(hi)
    while mean_prob(hi) < target_p:
        hi *= 2.0
        if hi > 1e12:  # extremely extreme; still keep going would be pointless
            break

    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        m = mean_prob(mid)
        if abs(m - target_p) <= tol:
            return mid
        if m < target_p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)





###### sparse graphs ######


# create sprase_graph
def create_sparse_graph(graph: nx.Graph, edges_length: Union[int, Dict[tuple, int]]):
    graph_copy = graph.copy()
    n = graph.number_of_nodes()
    if isinstance(edges_length, int):
        edges_length = {edge: edges_length for edge in graph.edges}
    for edge, l in edges_length.items():
        if l > 0:
            graph_copy.remove_edges_from([edge])
            new_edges = [(i, i+1) for i in range(n+1, n+l)]
            new_edges += [(edge[0], n+1), (edge[1], n+l)]
            graph_copy.add_edges_from(new_edges)
            n += l
    return graph_copy

# create sprase_graph with nodes positions equal position
def create_sparse_graph_with_pos(graph: nx.Graph, edges_length: dict):
    graph_copy = graph.copy()
    nodes_pos = nx.get_node_attributes(graph, "pos")
    if len(nodes_pos) != len(graph.nodes):
        raise ValueError("not all nodes have positions")
    for edge, l in edges_length.items():
        # Generate l equal nodes between edge[0], edge[1]
        u, v = np.array(nodes_pos[edge[0]]), np.array(nodes_pos[edge[1]])
        new_nodes = [
            k / (l + 1) * u + (1 - k / (l + 1)) * v
            for k in range(1, l + 1)
        ]        
        new_nodes = list(map(tuple, new_nodes))
        # Add the new edges
        if len(new_nodes) > 0:
            graph_copy.remove_edges_from([edge])
            new_edges = [(new_nodes[i], new_nodes[i + 1]) for i in range(l - 1)]
            new_edges += [(edge[1], new_nodes[0]), (edge[0], new_nodes[-1])]
            graph_copy.add_edges_from(new_edges)
    
    # Set positions
    pos = nx.get_node_attributes(graph, "pos")
    nx.set_node_attributes(
        graph_copy,
        {
            node: np.array(node) if node not in graph.nodes else pos[node]
            for node in graph_copy.nodes
        },
        "pos"
    )
    labels = {node: node if node in graph.nodes else i for i, node in enumerate(graph_copy.nodes)}
    graph_copy = nx.relabel_nodes(graph_copy, labels)
    return graph_copy


###### simple saidi ######


def all_subsets(x: list) -> list:
    """## Get a list and return a list of al its subsets

    ### Args:
        - `x (list)`: a list

    ### Returns:
        - `list`: a list of all x subsets
    """
    res = []
    for k in range(len(x) + 1):
        res.extend(combinations(x, k))
    return res

def path_graph_saidi(path_graph: nx.MultiGraph, source):
    # _validate_path_graph(path_graph, source)
    prob_class = next(iter(nx.get_edge_attributes(path_graph, "prob").values())).__class__
    # handle self edges
    deg1_nodes = [node for node, degree in path_graph.degree if degree == 1]
    if len(deg1_nodes) == 1:
        return prob_class([0.0])
    if path_graph.number_of_nodes() <= 2:
        return prob_class([0.0])
    one = prob_class([1.0])
    saidi = prob_class([0.0])
    node_failing_prob = prob_class([0.0])
    nodes_order = list(nx.dfs_preorder_nodes(path_graph, source=source))[:-1]
    prev_node = None
    for current_node in nodes_order:
        if prev_node is not None:
            edge_failing_prob = path_graph.edges[(prev_node, current_node, 0)]['prob']
            node_failing_prob = one - (one - node_failing_prob) * (one - edge_failing_prob)
            saidi += node_failing_prob * path_graph.nodes[current_node]['weight']
        prev_node = current_node
    total_weight = sum([path_graph.nodes[node]['weight'] for node in nodes_order[1:]])
    if total_weight == 0:
        print("Warning: total weight in path_Graph is zero")
        return prob_class([0.0])
    saidi /= total_weight
    return saidi

def path_graph_saidi_two_sources(path_graph: nx.MultiGraph, source):
    # _validate_path_graph(path_graph)
    path_graph_copy = path_graph.copy()
    if not nx.is_tree(path_graph_copy): # is cycle
        path_graph_copy = cycle_graph_into_path_with_two_sources(path_graph_copy, source)
    n = path_graph_copy.number_of_nodes()
    if n <= 2:
        return 0.0
    # define prob_class
    prob_class = next(iter(nx.get_edge_attributes(path_graph_copy, "prob").values())).__class__
    one = prob_class([1.0])
    saidi = prob_class([0.0])
    nodes_order = list(nx.dfs_preorder_nodes(path_graph_copy, source=source))
    prob_left = [prob_class([0]) for _ in range(n)]
    for i in range(1, n - 1):
        edge_failing_prob = path_graph_copy.edges[(nodes_order[i - 1], nodes_order[i], 0)]['prob']
        prob_left[i] = one - (one - prob_left[i - 1]) * (one - edge_failing_prob)
    prob_right = [prob_class([0]) for _ in range(n)]
    for i in range(n - 2, 0, -1):
        edge_failing_prob = path_graph_copy.edges[(nodes_order[i + 1], nodes_order[i], 0)]['prob']
        prob_right[i] = one - (one - prob_right[i + 1]) * (one - edge_failing_prob)
    weights = np.array([path_graph_copy.nodes[node]['weight'] for node in nodes_order])
    total_weight = weights[1:n-1].sum()
    if total_weight == 0:
        print("Warning: total weight in cycle graph is zero")
        print(weights)
        print(nx.get_node_attributes(path_graph_copy, 'weight'))
        return prob_class([0.0])
    saidi = np.sum(np.array(prob_left) * np.array(prob_right) * weights) / total_weight
    return saidi 

def cycle_graph_into_path_with_two_sources(path_graph: nx.MultiGraph, source: int) -> nx.MultiGraph:
    path_graph_copy = path_graph.copy()
    node1, node2, data = next(iter(path_graph_copy.edges(source, data=True)))
    new_node = max(path_graph.nodes) + 1
    path_graph_copy.remove_edge(node1, node2)
    new_node = max(path_graph.nodes) + 1
    path_graph_copy.add_edge(node2, new_node, **data)
    path_graph_copy.nodes[new_node]['weight'] = 0
    return path_graph_copy

def _validate_path_graph(path_graph: nx.MultiGraph, source=None):
    n = path_graph.number_of_nodes()
    degree1_nodes = [node for node, degree in path_graph.degree if degree == 1]
    degree2_nodes = [node for node, degree in path_graph.degree if degree == 2]
    if (len(degree1_nodes) != 2) or (len(degree2_nodes) != n - 2) or (not nx.is_connected(path_graph)):
        raise ValueError("graph is not a path graph")
    if source is not None and source not in degree1_nodes:
        raise ValueError(f"the source {source} is not in the ent pf the path graph")
