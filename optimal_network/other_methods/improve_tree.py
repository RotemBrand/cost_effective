import networkx as nx
import numpy as np
from itertools import combinations
from typing import List, Tuple
from tqdm import tqdm
from scipy.spatial.distance import pdist, squareform


import time
TIME = 0
def toc(s: str=''):
    global TIME
    if s != '':
        print(f'{s}: {time.time() - TIME}')
    TIME = time.time()

type EdgeList = List[Tuple]

def improve_tree(G_tree: nx.Graph, R: int, root, weight_quantile: float=1, show_progress: bool=False) -> Tuple[nx.Graph, EdgeList]:
    """
    Iteratively add chords to improve reliability by covering high-score edges.
    Weight = geometric distance between chord endpoints.
    Edge scores = number of nodes disconnected if edge fails.
    """
    chords = quantile_edges(G_tree, weight_quantile)
    # Compute initial edge scores (disconnection impact)
    edge_scores = _compute_edge_scores(G_tree, root)

    # Precompute chord coverage and weights
    coverage = _compute_chord_coverages(G_tree, root, chords)
    weight = _compute_chord_weights(G_tree, chords)

    selected = []
    scores = edge_scores.copy()

    progress = lambda x: tqdm(x, desc="improve tree") if show_progress else x
    for i in progress(list(range(R))):
        toc()
        best_chord, best_ratio = None, 0
        for c in chords:
            benefit = sum(scores[_key(e)] for e in coverage[_key(c)])
            if benefit <= 0:
                continue
            ratio = benefit / weight[_key(c)]
            if ratio > best_ratio:
                best_chord, best_ratio = c, ratio

        if best_chord is None:
            break

        selected.append(best_chord)
        # Remove covered edges (set their score to 0)
        for e in coverage[_key(best_chord)]:
            scores[_key(e)] = 0
    # return graph
    G_copy = G_tree.copy()
    G_copy.add_edges_from(selected)
    return G_copy, selected


def quantile_edges(G: nx.Graph, q: float) -> EdgeList:
    """Return node pairs whose Euclidean distance (from 'pos') is below the q-quantile. assume that the network has pos attribute"""
    pos = np.array(list(nx.get_node_attributes(G, "pos").values()))
    nodes = list(G.nodes)

    dist_mat = squareform(pdist(pos))

    threshold = np.quantile(dist_mat.flatten(), q)
    edges = [
        (nodes[i], nodes[j])
        for i, j in combinations(range(len(nodes)), 2)
        if dist_mat[i, j] < threshold
    ]
    return edges



def _key(e: tuple) -> tuple:
    return tuple(sorted(e))

# ------------------------------------------------------------
# 1. Compute subtree disconnection scores
# ------------------------------------------------------------
def _compute_edge_scores(G: nx.Graph, root):
    """
    Compute for each edge (u,v) in a tree the number of nodes
    disconnected if that edge fails.
    """
    edge_scores = {}
    n = len(G)

    def dfs(node, parent=None):
        size = 1
        for nei in G[node]:
            if nei == parent:
                continue
            sub_size = dfs(nei, node)
            edge_scores[_key((node, nei))] = sub_size
            size += sub_size
        return size

    dfs(root)
    return edge_scores


# ------------------------------------------------------------
# 2. Tree utilities
# ------------------------------------------------------------
def _compute_parent_depth(G: nx.Graph, root):
    """Compute parent and depth arrays for all nodes in a rooted tree."""
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


def _chord_cycle_edges(u, v, parent, depth):
    """Return the set of tree edges forming the unique cycle created by adding (u,v)."""
    uu, vv = u, v
    path_u, path_v = [uu], [vv]

    # move both nodes up until same depth
    while depth[uu] > depth[vv]:
        uu = parent[uu]
        path_u.append(uu)
    while depth[vv] > depth[uu]:
        vv = parent[vv]
        path_v.append(vv)

    # climb together until LCA found
    while uu != vv:
        uu = parent[uu]
        vv = parent[vv]
        path_u.append(uu)
        path_v.append(vv)

    # merge paths and convert to edges
    path = path_u + path_v[-2::-1]
    return {tuple(sorted((path[i], path[i+1]))) for i in range(len(path)-1)}


# ------------------------------------------------------------
# 3. Weight (distance) computation
# ------------------------------------------------------------
def _chord_length(G: nx.Graph, chord):
    """Compute Euclidean distance between chord endpoints using node 'pos'."""
    u, v = chord
    p1, p2 = np.array(G.nodes[u]["pos"]), np.array(G.nodes[v]["pos"])
    return np.linalg.norm(p1 - p2)


def _compute_chord_weights(G: nx.Graph, chords):
    """Return dict of weights (distance) for all candidate chords."""
    return {_key(c): _chord_length(G, c) for c in chords}


# ------------------------------------------------------------
# 4. Coverage computation
# ------------------------------------------------------------
def _compute_chord_coverages(G: nx.Graph, root, chords):
    """Precompute for each chord which tree edges it covers (its cycle)."""
    parent, depth = _compute_parent_depth(G, root)
    coverage = {_key(c): _chord_cycle_edges(*c, parent, depth) for c in chords}
    return coverage

