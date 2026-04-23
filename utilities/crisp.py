import numpy as np
import networkx as nx

def key(e: tuple) -> tuple:
    """Normalize edge tuple representation."""
    if len(e) == 2:
        return tuple(sorted(e))
    elif len(e) == 3:
        return (*sorted(e[:2]), e[2])
    else:
        raise ValueError(f"length of e should be 2 or 3, {e}")

def add_edges_to_random_tree(graph: nx.Graph, r: int, seed: int = 42, tree_factor: float = 1.0) -> nx.Graph:
    """
    Construct a robust network by generating a pseudo-random spanning tree
    and greedily adding r redundant chords to cover the maximum number of tree edges.
    """
    # get tree
    tree = random_spanning_tree(graph, seed=seed, factor=tree_factor)
    
    # get all cover paths
    chords = list(map(key, set(graph.edges) - set(tree.edges)))
    
    # precompute the cycle edges each chord creates
    root = next(iter(tree.nodes))
    parent, depth = _compute_parent_depth(tree, root)
    chords_cycles = {c: _chord_cycle_edges(c[0], c[1], parent, depth) for c in chords}
    
    # greedily select chords
    selected_chords = cover_tree_with_max_r_greedy(chords_cycles, n_edges=r)

    # build subgraph
    subgraph = tree.copy()
    subgraph.add_edges_from(selected_chords)
    
    # Use edge_subgraph to retain original edge attributes
    subgraph = nx.edge_subgraph(graph, subgraph.edges).copy()
    
    for e in subgraph.edges:
        subgraph.edges[e]["source"] = "tree" if e in tree.edges else "chord"

    return subgraph

def random_spanning_tree(graph: nx.Graph, seed: int = None, factor: float = 1.0) -> nx.Graph:
    """
    Generate a spanning tree with randomized edge weights.
    factor determines the amount of randomization (0 to 1).
    """
    if not (0 <= factor <= 1):
        raise ValueError(f"Factor should be 0 <= factor <= 1 instead of {factor}")
    
    graph_copy = graph.copy()
    rng = np.random.default_rng(seed)
    
    for e in graph_copy.edges:
        graph_copy.edges[e]['weight'] *= rng.random() * factor + (1 - factor)
        
    return nx.minimum_spanning_tree(graph_copy, weight="weight")

def _compute_parent_depth(G: nx.Graph, root: int) -> tuple[dict, dict]:
    """Compute parent and depth arrays for all nodes in a rooted tree using DFS/BFS."""
    parent = {root: None}
    depth = {root: 0}
    stack = [root]
    while stack:
        node = stack.pop()
        for nei in G[node]:
            if nei not in parent:
                parent[nei] = node
                depth[nei] = depth[node] + 1
                stack.append(nei)
    return parent, depth

def _chord_cycle_edges(u: int, v: int, parent: dict, depth: dict) -> set:
    """Return the set of tree edges forming the unique cycle created by adding (u,v)."""
    uu, vv = u, v
    edges = set()
    
    # Move both nodes up until they are at the same depth
    while depth[uu] > depth[vv]:
        p = parent[uu]
        edges.add(key((uu, p)))
        uu = p
    while depth[vv] > depth[uu]:
        p = parent[vv]
        edges.add(key((vv, p)))
        vv = p
        
    # Climb together until the Lowest Common Ancestor (LCA) is found
    while uu != vv:
        pu = parent[uu]
        pv = parent[vv]
        edges.add(key((uu, pu)))
        edges.add(key((vv, pv)))
        uu = pu
        vv = pv
        
    return edges

def cover_tree_with_max_r_greedy(chords_cycles: dict[tuple, set], n_edges: int) -> list[tuple]:
    """
    Greedily select n_edges chords that cover the maximum number of unique tree edges.
    """
    selected_chords = []
    # Copy sets to avoid mutating the input dictionary
    remaining_covers = {chord: set(edges) for chord, edges in chords_cycles.items()}
    
    for _ in range(n_edges):
        if not remaining_covers:
            break
            
        # Find the chord covering the most uncovered edges
        best_chord = max(remaining_covers, key=lambda c: len(remaining_covers[c]))
        best_cover = remaining_covers[best_chord]
        
        # If the best chord covers 0 new edges, we can stop early
        if len(best_cover) == 0:
            break
            
        selected_chords.append(best_chord)
        del remaining_covers[best_chord]
        
        # Remove newly covered edges from the remaining candidate chords
        for c in remaining_covers:
            remaining_covers[c] -= best_cover
            
    return selected_chords
