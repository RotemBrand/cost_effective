import numpy as np
import networkx as nx
import pandas as pd
from typing import Tuple, List, Union
from tqdm import tqdm

# Define a custom type hint for an edge representing a tuple of integers
Edge = Tuple[int, ...]

def key(e: Edge) -> Edge:
    """
    Standardize the representation of an edge so that undirected edges
    always have the same tuple representation, regardless of node order.

    Args:
        e: A tuple representing an edge. Usually (u, v) or (u, v, key)

    Returns:
        A sorted tuple representing the edge.
    """
    if len(e) == 2:
        return tuple(sorted(e))
    elif len(e) == 3:
        return (*sorted(e[:2]), e[2])
    else:
        raise ValueError(f"Length of e should be 2 or 3, but got {e}")


def random_spanning_tree(graph: nx.Graph, seed: Union[int, None] = None, factor: float = 1.0) -> nx.Graph:
    """
    Generate a random spanning tree by slightly perturbing the edge weights.

    Args:
        graph: The base graph to extract a spanning tree from.
        seed: Random seed for reproducibility.
        factor: A float between 0 and 1 controlling the randomness. A factor of 0
                yields the exact minimum spanning tree, while 1 yields a highly randomized tree.

    Returns:
        A minimum spanning tree of the graph with perturbed weights.
    """
    if factor < 0 or factor > 1:
        raise ValueError(f"Factor should be 0 <= factor <= 1 instead of {factor}")
    
    graph_copy = graph.copy()
    rng = np.random.default_rng(seed)
    
    # Perturb the weight of each edge based on the factor
    for e in graph_copy.edges:
        graph_copy.edges[e]['weight'] *= (rng.random() * factor) + (1 - factor)
        
    return nx.minimum_spanning_tree(graph_copy, weight="weight")


def all_paths_tree(tree: nx.Graph, st_list: list) -> List[List[int]]:
    """
    Find the unique path within a tree between pairs of source-target nodes.
    
    This function uses BFS to assign a hierarchical "sequence" (like a Dewey Decimal index)
    to each node. To find the path between two nodes, it finds their lowest common ancestor
    by comparing their sequences, and reconstructs the path upwards and then downwards.

    Args:
        tree: The spanning tree graph.
        st_list: A list of source-target tuples (e.g. chords) to find paths for.

    Returns:
        A list of paths, where each path is a list of node IDs.
    """
    # 1. Assign a sequence to each node using BFS
    root = next(iter(tree.nodes))
    tree.nodes[root]['seq'] = tuple([0])
    
    for node, neighs in nx.bfs_successors(tree, root):
        node_seq = tree.nodes[node]['seq']
        for i, neigh in enumerate(neighs):
            tree.nodes[neigh]["seq"] = tuple(list(node_seq) + [i])
            
    # Reverse mapping from sequence to node
    seq_to_node = {
        tuple(seq): node for node, seq in nx.get_node_attributes(tree, "seq").items()
    }

    # 2. Reconstruct paths using the sequences
    paths = []
    for s, t in tqdm(st_list, desc="Calculating tree paths"):
        seq1, seq2 = tree.nodes[s]["seq"], tree.nodes[t]["seq"]
        
        # Find the index where the sequences diverge (Lowest Common Ancestor)
        i = 0
        while i < len(seq1) and i < len(seq2) and seq1[i] == seq2[i]:
            i += 1

        # Path goes up from source to LCA, then down from LCA to target
        path = [seq_to_node[seq1[:j]] for j in range(len(seq1), i, -1)]
        path += [seq_to_node[seq2[:j]] for j in range(i, len(seq2) + 1)]
        paths.append(path)
        
    return paths


def cover_tree_with_max_r_gridy(graph: nx.Graph, chords_cycle: dict, tree: nx.Graph, n_edges: int) -> List[Edge]:
    """
    Greedily select a given number of chords to maximize the number of tree edges 
    covered by their resulting fundamental cycles.

    Args:
        graph: The base graph.
        chords_cycle: A dictionary mapping a chord edge to the list of nodes forming its fundamental cycle.
        tree: The spanning tree.
        n_edges: The number of chords to select (r).

    Returns:
        A list of selected chord edges.
    """
    selected_chords = []
    
    # Create a DataFrame where each chord has a set of the nodes/edges it covers
    chord_cycle_data = pd.DataFrame(
        [[k, set(v)] for k, v in chords_cycle.items()], 
        columns=["chord", "cover_edges"]
    )
    
    # Greedily pick the chord that covers the maximum number of new elements
    for _ in tqdm(list(range(n_edges)), desc="Greedily adding chords"):
        # Calculate how many elements each chord covers
        chord_cycle_data["cover_edges_len"] = chord_cycle_data["cover_edges"].apply(len)
        
        # Select the chord that covers the most elements
        best_row_idx = chord_cycle_data["cover_edges_len"].idxmax()
        chord, cover_edges, _ = chord_cycle_data.loc[best_row_idx]
        
        # Remove the newly covered elements from all remaining chords
        chord_cycle_data["cover_edges"] = chord_cycle_data["cover_edges"].apply(lambda x: x - cover_edges)
        selected_chords.append(chord)
    
    return selected_chords


def add_edges_to_random_tree(graph: nx.Graph, r: int, seed: int = 42, tree_factor: float = 1.0) -> nx.Graph:
    """
    Generate a robust network design by starting with a randomized spanning tree and
    greedily adding 'r' chords that maximize structural coverage.

    Args:
        graph: The base graph containing all possible edges.
        r: The number of redundant edges (chords) to add.
        seed: Random seed for reproducibility.
        tree_factor: Controls how randomized the initial spanning tree is.

    Returns:
        A robust subgraph containing the tree and the selected chords.
    """
    # 1. Generate a randomized base spanning tree
    tree = random_spanning_tree(graph, seed=seed, factor=tree_factor)
    
    # 2. Identify all possible chords (edges not in the tree)
    chords = list(map(key, set(graph.edges) - set(tree.edges)))
    
    # 3. For each chord, find its fundamental cycle (the unique path in the tree)
    chords_cycles = dict(zip(chords, all_paths_tree(tree, chords)))
    
    # 4. Greedily select 'r' chords that cover the maximum amount of the tree
    selected_chords = cover_tree_with_max_r_gridy(graph, chords_cycles, tree, n_edges=r)

    # 5. Build the final subgraph
    subgraph = tree.copy()
    subgraph.add_edges_from(selected_chords)
    subgraph = nx.edge_subgraph(graph, subgraph.edges)
    
    # Tag edges for visualization or analysis
    for e in subgraph.edges:
        subgraph.edges[e]["source"] = "tree" if e in tree.edges else "chord"

    return subgraph
