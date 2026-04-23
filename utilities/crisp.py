import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Literal,  Hashable, Optional, Set
import gurobipy as gp
from gurobipy import GRB
from utilities.helper import draw_network
from indexes.utilities import get_skeleton_graph
from icecream import ic
import heapq
from collections import deque
from tqdm import tqdm
import time
import pandas as pd

TIME = 0
def toc(s: str=''):
    global TIME
    if s != '':
        print(f'{s}: {time.time() - TIME}')
    TIME = time.time()


Edge = Tuple[int]

def key(e):
    if len(e) == 2:
        return tuple(sorted(e))
    elif len(e) == 3:
        return (*sorted(e[:2]), e[2])
    else:
        raise ValueError(f"length of e should be 2 or 3, {e}")

def crisp2_connected(graph: nx.Graph, chain_min_len: int=0, max_r: int=None) -> nx.Graph:
    # mst
    mst: nx.Graph = prim_mst(graph, chain_min_len=chain_min_len)
    one_deg_nodes = [n for n, d in mst.degree if d == 1]
    all_paths = dict(nx.all_pairs_dijkstra_path(mst, weight="weight"))
    chords_cycles = {}
    for v, data_v in all_paths.items():
        for u, path in data_v.items():
            if v <= u and not mst.has_edge(v, u) and graph.has_edge(v, u):
                chords_cycles[key((int(u), int(v)))] = list(map(int, path))
    print(list(chords_cycles.items())[:2])
    # for e in graph.edges:
    #     if mst.has_edge(*e):
    #         continue
    #     cycle = all_paths[e[0]][e[1]]
    #     chords_cycles[key(e)] = cycle
    
    # ILP
    selected_chords = cover_ilp(graph, chords_cycles, mst, max_edges=max_r)
    subgraph = mst.copy()
    subgraph.add_edges_from(selected_chords)
    subgraph = nx.edge_subgraph(graph, subgraph.edges)

    return subgraph



def get_edge_to_cover_chords(chords_cycle: dict[Edge, List]) -> Dict[Tuple, set]:
    edge_to_cover: Dict[Tuple, set] = {}
    for chord, cycle in chords_cycle.items():
        for e in zip(cycle, np.roll(cycle, 1)):
            e = key(e)
            if e not in edge_to_cover:
                edge_to_cover[e] = set()
            edge_to_cover[e].add(key(chord))
    return edge_to_cover

def get_node_to_cover_chords(chords_cycle: dict[Edge, List]) -> Dict[Tuple, set]:
    node_to_cover: Dict[Tuple, set] = {}
    for chord, cycle in chords_cycle.items():
        for node in cycle:
            if node not in node_to_cover:
                node_to_cover[node] = set()
            node_to_cover[node].add(key(chord))
    return node_to_cover

def cover_ilp(graph: nx.Graph, chords_cycle: dict[Edge, List], mst: nx.Graph, max_edges: int=None):
    m = gp.Model("node_cover_by_cycles")
    m.Params.OutputFlag = 0
    
    x = m.addVars(list(chords_cycle.keys()), vtype=GRB.BINARY, name="x")

    # Cover every edge by at least one selected cycle
    edge_to_cover = get_edge_to_cover_chords(chords_cycle)
    node_to_cover = get_node_to_cover_chords(chords_cycle)
    if max_edges:
        m.addConstr(gp.quicksum(x[e] for e in chords_cycle) <= max_edges, name=f"sum")
    for edge in mst.edges:
        chords = edge_to_cover[key(edge)]
        m.addConstr(gp.quicksum(x[chord] for chord in chords) >= 1, name=f"edgecov[{edge}]")

    for node, chords in node_to_cover.items():
        m.addConstr(gp.quicksum(x[chord] for chord in chords) >= 1, name=f"nodecov[{node}]")

    # Minimize total chord weight
    m.setObjective(gp.quicksum(graph.edges[chord]["weight"] * x[chord] for chord in chords_cycle), GRB.MINIMIZE)
    m.optimize()

    status = m.Status
    if status != GRB.Status.OPTIMAL:
        raise ValueError("dont found optimal solution")

    selected = [e for e in chords_cycle if x[e].X > 0.5]
    return selected



def prim_mst(
    G: nx.Graph,
    start: Optional[Hashable] = None,
    weight: str = "weight",
    chain_min_len: int = 0,
) -> nx.Graph:
    """
    Prim's MST for a connected, undirected simple graph with a weight on every edge,
    with an additional spacing rule between fork nodes.

    Parameters
    ----------
    G : nx.Graph
        Connected, undirected, simple graph. Every edge must have the 'weight' attribute.
    start : optional node
        Starting node for Prim's growth. If None, an arbitrary node is chosen.
    weight : str
        Edge attribute name for weights.
    chain_min_len : int
        Minimum allowed hop-distance between **distinct fork nodes** (nodes whose degree in the
        growing tree reaches >= 3 for the first time). If 0, no restriction.

    Strategy for chain_min_len
    --------------------------
    Maintain:
      - `dist_fork`: minimal hop-distance from each node to the nearest *current* fork in H
        (∞ if none yet). Updated incrementally.
      - `forbidden`: nodes with `dist_fork[node] <= chain_min_len`. Nodes in this set may be added
        to the tree as chains, but **may not become new forks**.

    When we add a new fork f, we run a BFS *on the current tree H* up to `chain_min_len` hops to
    update `dist_fork` and `forbidden`. When we add a new node v via edge (u,v), we set
    `dist_fork[v] = min(dist_fork[v], dist_fork[u] + 1)` and, if within radius, add v to
    `forbidden`. This ensures future fork checks see v as forbidden if it's too close to an
    existing fork, even if v was added after that fork was created.

    Returns
    -------
    H : nx.Graph
        The resulting minimum spanning tree that also respects the fork-spacing constraint.
    """
    if G.is_directed():
        raise ValueError("Prim requires an undirected graph.")
    if G.number_of_nodes() == 0:
        return nx.Graph()

    if start is None:
        start = next(iter(G.nodes))

    H = nx.Graph()
    H.add_nodes_from(G.nodes(data=True))

    in_tree = {start}
    forks: Set[Hashable] = set()
    forbidden: Set[Hashable] = set()  # nodes too close (<= chain_min_len hops) to any fork
    dist_fork: Dict[Hashable, int] = {n: 10**9 for n in H.nodes}  # ∞ as large int

    def _propagate_from_fork(f: Hashable):
        """BFS from newly created fork f over current H up to chain_min_len hops;
        update dist_fork and forbidden incrementally."""
        if chain_min_len <= 0:
            return
        dq = deque([f])
        if dist_fork[f] > 0:
            dist_fork[f] = 0
        # ensure the fork itself is forbidden
        forbidden.add(f)
        while dq:
            u = dq.popleft()
            du = dist_fork[u]
            if du >= chain_min_len:
                continue
            for w in H.neighbors(u):
                nd = du + 1
                if nd < dist_fork[w]:
                    dist_fork[w] = nd
                    if nd <= chain_min_len:
                        forbidden.add(w)
                    dq.append(w)

    # priority queue of candidate edges (weight, u_in_tree, v_outside)
    pq = []
    for v, data in G[start].items():
        heapq.heappush(pq, (data[weight], start, v))

    while pq and len(in_tree) < G.number_of_nodes():
        w, u, v = heapq.heappop(pq)
        if v in in_tree:
            continue
        # Check spacing rule only if adding (u,v) would create a *new* fork at u
        if chain_min_len > 0:
            deg_u = H.degree(u)
            will_make_new_fork = (deg_u >= 2) and (u not in forks)
            if will_make_new_fork and dist_fork[u] < chain_min_len:
                # skip this edge; try another candidate
                continue

        # Accept edge u-v
        H.add_edge(u, v, **dict(G[u][v]))  # copy all edge attrs (incl. weight)
        in_tree.add(v)

        # Incremental update for the new node's distance-from-fork
        if chain_min_len > 0:
            if dist_fork[u] + 1 < dist_fork[v]:
                dist_fork[v] = dist_fork[u] + 1
            if dist_fork[v] <= chain_min_len:
                forbidden.add(v)

        # If this made u a *new* fork (deg became 3), register and propagate
        if chain_min_len > 0 and H.degree(u) >= 3 and u not in forks:
            forks.add(u)
            _propagate_from_fork(u)

        # Push new frontier edges from v
        for x, d in G[v].items():
            if x not in in_tree:
                heapq.heappush(pq, (d[weight], v, x))

    # sanity: MST must have |V|-1 edges if G is connected
    if H.number_of_edges() != G.number_of_nodes() - 1:
        raise RuntimeError(
            "Could not span all nodes under the fork-spacing constraint; graph or constraint may be infeasible."
        )
    return H


###### ILP ######
def add_edges_to_random_tree(graph: nx.Graph, r: int, seed: int=42, tree_factor: int=1) -> nx.Graph:
    # get tree
    tree = random_spanning_tree(graph, seed=seed, factor=tree_factor)
    # get all cover paths
    chords = list(map(key, set(graph.edges) - set(tree.edges)))
    chords_cycles = dict(zip(chords, all_paths_tree(tree, chords)))
    # ilp
    selected_chords = cover_tree_with_max_r_gridy(graph, chords_cycles, tree, n_edges=r)

    # return subgraph
    subgraph = tree.copy()
    subgraph.add_edges_from(selected_chords)
    subgraph = nx.edge_subgraph(graph, subgraph.edges)
    for e in subgraph.edges:
        subgraph.edges[e]["source"] = "tree" if e in tree.edges else "chord"

    return subgraph

def random_spanning_tree(graph: nx.Graph, seed: int=None, factor: int=1):
    if factor < 0 or factor > 1:
        raise ValueError(f"Factor should be 0 < factor < 1 insted of {factor}")
    graph_copy = graph.copy()
    rng = np.random.default_rng(seed)
    for e in graph_copy.edges:
        graph_copy.edges[e]['weight'] *= rng.random() * factor + (1 - factor)
    return nx.minimum_spanning_tree(graph_copy, "weight")

def cover_tree_with_max_r(graph: nx.Graph, chords_cycle: dict[Edge, List], tree: nx.Graph, n_edges: int):
    m = gp.Model("node_cover_by_cycles")
    m.Params.OutputFlag = 0
    
    edges = list(map(key, tree.edges))
    x = m.addVars(list(chords_cycle.keys()), vtype=GRB.BINARY, name="x")
    E = m.addVars(edges, vtype=GRB.BINARY, name="E")

    # Cover every edge by at least one selected cycle
    edge_to_cover = get_edge_to_cover_chords(chords_cycle)
    m.addConstr(gp.quicksum(x[e] for e in chords_cycle) == n_edges, name=f"sum")
    for e in edges:
        chords = edge_to_cover[key(e)]
        m.addConstr(gp.quicksum(x[chord] for chord in chords) >= E[e], name=f"edgecov[{e}]")


    # Minimize total chord weight
    m.setObjective(gp.quicksum(E[e] for e in edges), GRB.MAXIMIZE)
    m.optimize()

    status = m.Status
    if status != GRB.Status.OPTIMAL:
        raise ValueError("dont found optimal solution")

    selected = [e for e in chords_cycle if x[e].X > 0.5]
    return selected

def cover_tree_with_max_r_gridy(graph: nx.Graph, chords_cycle: dict[Edge, List], tree: nx.Graph, n_edges: int):
    # edge_to_covers = get_edge_to_cover_chords(chords_cycle)

    selected_chords = []
    chord_cycle_data = pd.DataFrame([[k, set(v)] for k, v in chords_cycle.items()], columns=["chord", "cover_edges"])
    for i in tqdm(list(range(n_edges)), desc="add chords"):
        chord_cycle_data["cover_edges_len"] = chord_cycle_data["cover_edges"].apply(len)
        chord, cover_edges, _ = chord_cycle_data.loc[chord_cycle_data["cover_edges_len"].idxmax()]
        chord_cycle_data["cover_edges"] = chord_cycle_data["cover_edges"].apply(lambda x: x - cover_edges)
        selected_chords.append(chord)
    

    return selected_chords

def all_paths_tree(tree: nx.Graph, st_list: list):
    # give seq for each node
    root = next(iter(tree.nodes))
    tree.nodes[root]['seq'] = tuple([0])
    for node, neighs in nx.bfs_successors(tree, root):
        node_seq = tree.nodes[node]['seq']
        for i, neigh in enumerate(neighs):
            tree.nodes[neigh]["seq"] = tuple(list(node_seq) + [i])
    seq_to_node = {
        tuple(seq): node for node, seq in nx.get_node_attributes(tree, "seq").items()
    }

    # paths
    paths = []
    for s, t in tqdm(st_list, desc="paths"):
        seq1, seq2 = tree.nodes[s]["seq"], tree.nodes[t]["seq"]
        # meet
        i = 0
        while i < len(seq1) and i < len(seq2) and seq1[i] == seq2[i]:
            i += 1

        path = [seq_to_node[seq1[:j]] for j in range(len(seq1), i, -1)]
        path += [seq_to_node[seq2[:j]] for j in range(i, len(seq2) + 1)]
        paths.append(path)
    return paths


###### ILP ######

def key_edges(edges: Tuple[Edge]):
    return tuple(sorted([key(e) for e in edges]))

def total_score(graph: nx.Graph) -> float:
    return sum(nx.get_edge_attributes(graph, "weight").values())

def get_all_chains(graph: nx.Graph):
    sources = [node for node, d in graph.degree if d >= 3]
    if len(sources) == 0:
        sources = [next(iter(graph.nodes))]
    chains = _get_all_chains_from_graphs(graph, source=next(iter(sources)))
    return chains

def _get_all_chains_from_graphs(graph: nx.Graph):
    sources = [node for node, d in graph.degree if d >= 3]
    if len(sources) == 0:
        sources = [next(iter(graph.nodes))]
    chains = _get_all_chains_with_source(graph, source=next(iter(sources)))
    return chains

def _get_all_chains_with_source(subgraph: nx.Graph, source: int) -> List[List[int]]:
    gi = GraphImprover(subgraph, sources=[source])
    chains = [
        gi.chain(chain).nodes[1:-1]
        for chain in gi.structure.edges
    ]
    # create a graph of chains in the same mcs
    chains_graph = nx.Graph()
    for mcs in gi.mcs_list:
        chains_graph.add_edge(*[chain.name for chain in mcs.chains])
    # draw_network(chains_graph)
    # plt.show()
    # add cliques of the chain graph to the chains list
    cliques = list(nx.find_cliques(chains_graph))
    nodes_order = nx.dfs_preorder_nodes(subgraph, source)
    for clique in cliques:
        nodes_set = set([
            node for chain_name in clique
            for node in gi.chain(chain_name).nodes[1:-1]
        ])
        nodes = [node for node in nodes_order if node in nodes_set]
        chains.append(nodes)
        chains_graph.remove_edges_from(list(combinations(clique, 2)))
    # add other mcses to the chains list
    for chain_tuple in chains_graph.edges:
        nodes = [
            node for chain_name in chain_tuple
            for node in gi.chain(chain_name).nodes[1:-1]
        ]
        chains.append(nodes)
    return chains



def get_max_chain_len(graph: nx.Graph) -> float:
    if len(list(nx.k_edge_components(graph, 2))) == 1:
        chains = get_all_chains(graph)
        return max(map(len, chains))
    return np.inf

def draw_2_comps(graph: nx.Graph, **args):
    comps = list(filter(lambda x: len(x) > 1, nx.k_edge_components(graph, 2)))
    import seaborn as  sns


    text = f"nodes: {graph.number_of_nodes()}, n_comps = {len(comps)}, weight = {total_score(graph) / len(graph.nodes):.4f}, r = {len(graph.edges) - len(graph.nodes) + 1}"
    if len(comps) == 1:
        chains = get_all_chains(graph)
        max_chain = max([len(chain) for chain in chains])
        text += f", max_chain: {max_chain}"
        colors = sns.color_palette(n_colors=len(chains))
        node_to_color = {node: color for chain, color in zip(chains, colors) for node in chain}
    else:
        colors = sns.color_palette(n_colors=len(comps))
        node_to_color = {node: color for comp, color in zip(comps, colors) for node in comp}
    plt.gca().set_title(text)
    draw_network(
        graph,
        node_color=[node_to_color.get(node, "black") for node in graph.nodes],
        **args
    )


def minimal_2_degree_net(graph: nx.Graph, max_r: int=None, max_d=None, n_init_iters: int=3, n_max_chain_iters: int=3, plot_debug: bool=False) -> nx.Graph:
    # run first model with no edge favor
    m, edges, x, constrians = _init_net2_model(graph=graph, max_r=max_r, plot_debug=plot_debug, edges_to_favor=[])
    for it in range(n_init_iters):
        is_2connected, max_chain_len, subgraph = _add_comps_constrains_and_run(m, graph=graph, edges=edges, x=x, constrians=constrians, favor_existing_edges=False, max_d=max_d)
        if is_2connected:
            break
        if plot_debug:
            plt.figure()
            draw_2_comps(subgraph, chain_size=7, strc_size=30, with_labels=False)
            plt.show()

    if is_2connected:
        return subgraph

    # run second model with  edge favor
    best_subgraph, best_max_chain_len = None, np.inf
    selected = [e for e in edges if x[e].X > 0.5]
    m2, edges, x2, constrians2 = _init_net2_model(graph=graph, max_r=max_r, plot_debug=plot_debug, edges_to_favor=selected)
    n_examples_of_2_connected = 0
    while not is_2connected or (max_chain_len >= max_d):
        is_2connected, max_chain_len, subgraph = _add_comps_constrains_and_run(m2, graph=graph, edges=edges, x=x2, constrians=constrians2, favor_existing_edges=True, max_d=max_d)
        if is_2connected:
            n_examples_of_2_connected += 1
            print(f"found 2connected graph with max chain {max_chain_len}, it={n_examples_of_2_connected}, optimal len = {best_max_chain_len}")
            if max_chain_len < best_max_chain_len:
                best_max_chain_len = max_chain_len
                best_subgraph = subgraph
            if n_examples_of_2_connected >= n_max_chain_iters:
                break
        if plot_debug:
            plt.figure()
            draw_2_comps(subgraph, chain_size=7, strc_size=30, with_labels=False)
            plt.show()
    
    # return graph
    return best_subgraph


def _add_comps_constrains_and_run(m, graph: nx.Graph, edges: List[Edge], x, constrians: set, favor_existing_edges: bool, max_d: int):
    # find comps
    selected = [e for e in edges if x[e].X > 0.5]
    subgraph = nx.edge_subgraph(graph, selected).copy()
    comps = list(nx.k_edge_components(subgraph, 2))
    is_2connected = len(comps) == 1
    # add comps constrains
    if not is_2connected:
        for comp in comps:
            if len(comp) == 1:
                continue
            out_edges = key_edges([e for e in edges if len(set(e) & set(comp)) == 1])
            if out_edges not in constrians:
                m.addConstr(gp.quicksum(x[e] for e in out_edges) >= 2)
                constrians.add(out_edges)

    # add chains_constrain
    elif max_d is not None:
        sources = [node for node, d in subgraph.degree if d >= 3]
        if len(sources) == 0:
            sources = [next(iter(subgraph.nodes))]
        chains = ILP.get_all_chains(subgraph, source=next(iter(sources)))
        nodes_in_risk = set()
        for chain in chains:
            if len(chain) > max_d:
                # ic(str(chain))
                for i in range(len(chain) - max_d):
                    nodes = set(chain[i:i + max_d])
                    nodes_in_risk |= nodes
                    n_comps = nx.number_connected_components(nx.subgraph(subgraph, nodes))
                    subgraph_no_chains = subgraph.copy()
                    subgraph_no_chains.remove_edges_from([e for e in edges if set(e) & nodes])
                    node_to_comp = {
                        node: i
                        for i, comp in enumerate(nx.connected_components(subgraph_no_chains))
                        for node in comp 
                    }
                    # ic(len(nodes), sum([x[e].X for e in E if set(e) & nodes]))
                    m.addConstr(gp.quicksum(x[e] for e in edges if node_to_comp[e[0]] != node_to_comp[e[1]]) >= max_d + n_comps + 1)
        max_chain_len = get_max_chain_len(subgraph)
        if max_chain_len <= max_d:
            return True, max_chain_len, subgraph
    
    else:
        return True, max_chain_len, subgraph

    # define the new objective
    if favor_existing_edges:
        weights = {e: graph.edges[e]["weight"] for e in edges}
        for e in selected:
            weights[e] /= 100
        m.setObjective(gp.quicksum(weights[e] * x[e] for e in edges), GRB.MINIMIZE)

    # run optimization
    m.optimize()
    selected = [e for e in edges if x[e].X > 0.5]
    subgraph = nx.edge_subgraph(graph, selected).copy()
    is_2connected = len(list(nx.k_edge_components(subgraph, 2))) == 1
    max_chain_len = get_max_chain_len(subgraph)
    return is_2connected, max_chain_len, subgraph

def _init_net2_model(graph: nx.Graph, max_r: int=None, edges_to_favor: List[Edge]=[], plot_debug: bool=False):
    m = gp.Model("node_cover_by_cycles")
    m.Params.OutputFlag = 0
    
    edges = list(map(key, graph.edges))
    weights = {e: graph.edges[e]["weight"] for e in graph.edges}
    for e in edges_to_favor:
        weights[e] /= 100
    constrians = set()
    x = m.addVars(edges, vtype=GRB.BINARY, name="x")

    # node constrain 
    for node in graph.nodes:
        node_edges = key_edges(map(key, graph.edges(node)))
        constrians.add(node_edges)
        m.addConstr(gp.quicksum(x[e] for e in node_edges) >= 2, name=f"node[{node}]")

    # r_constrain
    if max_r is not None:
        max_num_edges = graph.number_of_nodes() + max_r - 1
        m.addConstr(gp.quicksum(x[e] for e in edges) <= max_num_edges, name="max_r")
    # Minimize total chord weight
    m.setObjective(gp.quicksum(weights[e] * x[e] for e in edges), GRB.MINIMIZE)
    m.optimize()

    status = m.Status
    if status != GRB.Status.OPTIMAL:
        return None
    return m, edges, x, constrians

def merge_chains(graph: nx.Graph, subgraph: nx.Graph) -> nx.Graph:
    if not is_graph_2connected(subgraph):
        return subgraph.copy()
    
    # extract all 4 cycles in the graph
    all_4_cycles = get_all_4_cycles(graph)
    all_chains = get_all_chains(subgraph)
    edge_to_chain = {
        key(edge): i
        for i, chain in enumerate(all_chains)
        for edge in subgraph.subgraph(chain).edges
    }
    node_order_in_chain = {
        node: (i, len(chain) - i)
        for chain in all_chains
        for i, node in enumerate(chain)
    }

    # filter cycles that have two oposite sides in different chains and 2 oposite edges wth no edge
    relevent_cycles = []
    for cycle in all_4_cycles:
        cycle_edges_chains = [edge_to_chain.get(key(e), -1) for e in cycle]
        is_exist_in_subgraph = [subgraph.has_edge(*e) for e in cycle]
        if (cycle_edges_chains[0] != cycle_edges_chains[2]) and ((is_exist_in_subgraph[1] == is_exist_in_subgraph[3] == -1)):
            relevent_cycles.append(cycle)
        elif (cycle_edges_chains[1] != cycle_edges_chains[3]) and ((is_exist_in_subgraph[0] == is_exist_in_subgraph[2] == -1)):
            relevent_cycles.append(np.roll(cycle, 1))
    
    # calculate the score for each merging
    weight_f = lambda e: graph.edges[e]["weight"]
    for cycle in relevent_cycles:
        weight_diff = weight_f(cycle[1]) + weight_f(cycle[3]) - weight_f(cycle[0]) - weight_f(cycle[1])
        chain1, chain2 = edge_to_chain[key(cycle[0])], edge_to_chain[key(cycle[2])]
        node_order1, node_order2 = node_order_in_chain[cycle[0][0]], node_order_in_chain[cycle[3][1]]


def is_contain_edges(graph: nx.Graph, edges: List[Edge], contains: List[bool]):
    """ Return True if for each edge in edges, the contains entry represent is the edge in the graph"""
    if len(contains) != len(edges):
        raise ValueError("contains and edges should have the same len")
    for e, contain in zip(edges, contains):
        if graph.has_edge(*e) != contain:
            return False
    return True

def redundancy(graph: nx.Graph):
    return graph.number_of_edges() - graph.number_of_nodes() + 1

def is_graph_2connected(graph: nx.Graph):
    return len(list(nx.k_edge_components(graph, 2))) == 1

def merge_comps(graph: nx.Graph, subgraph: nx.Graph, max_r: int):
    subgraph_copy = subgraph.copy()
    # build comps_graph
    comps = list(filter(
        lambda comp: len(comp) > 1,
        nx.k_edge_components(subgraph, 2)
    ))
    print(f"start merging {len(comps)} components")
    if len(comps) == 1:
        return subgraph
    comps_graph = nx.MultiGraph()
    for i, comp in enumerate(comps):
        comps_graph.add_node(i, comp=comp, size=len(comp))
    edge_to_comp = {
        key(edge): i
        for i, comp in enumerate(comps)
        for edge in subgraph.subgraph(comp).edges
    }
    j = min(enumerate(comps), key=lambda x: len(x[1]))[0]
    print({e for e, i in edge_to_comp.items() if i == j}, comps[j])
    node_to_comp = {
        node: i
        for i, comp in enumerate(comps)
        for node in comp
    }
    # check for merging options
    cycles4 = get_all_4_cycles(graph)
    for edges in cycles4:
        edge_comps = [edge_to_comp.get(key(e), -1) for e in edges]
        e0, e1, e2, e3 = tuple(map(key, edges))
        if edge_comps[0] != edge_comps[2] and edge_comps[1] == -1 and edge_comps[3] == -1 and edge_comps[0] != -1 and edge_comps[2] != -1:
            # ic(edge_comps[0], edge_comps[2], (e0, e2), (e1, e3))
            comps_graph.add_edge(edge_comps[0], edge_comps[2], old_edges=(e0, e2), new_edges=(e1, e3))
        if edge_comps[1] != edge_comps[3] and edge_comps[0] == -1 and edge_comps[2] == -1 and edge_comps[1] != -1 and edge_comps[3] != -1:
            # ic(edge_comps[1], edge_comps[3], (e1, e3), (e0, e2))
            comps_graph.add_edge(edge_comps[1], edge_comps[3], old_edges=(e1, e3), new_edges=(e0, e2))
    for e, d in comps_graph.edges.items():
        comps_graph.edges[e]["weight"] = (
            graph.edges[d["new_edges"][0]]["weight"] + graph.edges[d["new_edges"][1]]["weight"] - 
            graph.edges[d["old_edges"][0]]["weight"] + graph.edges[d["old_edges"][1]]["weight"]
        )
        # comp_size = [comps_graph.nodes[a]["size"] for a in (u, v)]
        # comps_graph.edges[e]["risk_diff"] = (
        #     (comp_size[0] + comp_size[1]) ** 2 - comp_size[0] ** 2 - comp_size[1] ** 2
        # )
        # comps_graph.edges[e]["score"] = comps_graph.edges[e]["risk_diff"] / comps_graph.edges[e]["weight"]  

    # # check for adding edge option
    # if redundancy(subgraph) < max_r:
    #     for e in graph.edges:
    #         comp1, comp2 = node_to_comp.get(e[0], -1), node_to_comp.get(e[1], -1) 
    #         if comp1 != comp2
    #             comps_graph.add_edge(comp1, comp2, weight=graph.edges[e]["weight"], )
            


    # start merging
    it = 0
    while comps_graph.number_of_nodes() != 1:
        it += 1
        print(comps_graph, it)
        mergings = comps_graph.edges.items()
        (comp1, comp2, _), data = max(mergings, key=lambda x: x[1]["weight"])
        comp1, comp2 = sorted([comp1, comp2])
        subgraph_copy.add_edges_from(data["new_edges"])
        subgraph_copy.remove_edges_from(data["old_edges"])
        # comps_graph.nodes[comp1]["comp"].extend(comps_graph.nodes[comp2]["comp"])
        # comps_graph.nodes[comp1]["size"] += comps_graph.nodes[comp2]["size"]
        print(f"merge {comp1} and {comp2}")
        nx.contracted_nodes(comps_graph, comp1, comp2, self_loops=False, copy=False)
        plt.figure()
        edge_color = ['red' if key(e) in data["new_edges"] else 'grey' for e in subgraph_copy.edges]
        draw_2_comps(
            subgraph_copy, chain_size=20, strc_size=60, with_labels=False, edge_color=edge_color
        )
        plt.show()
        if it > 300:
            break
    return subgraph_copy


        


from itertools import combinations
def get_all_4_cycles(graph: nx.Graph) -> List[List[Edge]]:
    graph_copy: nx.Graph = graph.copy()
    cycles = []
    for v in graph.nodes:
        neighbors = graph_copy[v]
        neighbors2 = {neigh: set(graph_copy[neigh]) for neigh in neighbors}
        for n1, n2 in combinations(neighbors, 2):
            common = neighbors2[n1] & neighbors2[n2] - set([v, n1, n2])
            for n3 in common:
                cycle = (v, n1, n3, n2)
                cycle_edges = list(zip(cycle, np.roll(cycle, 1)))
                cycles.append(cycle_edges)
        graph_copy.remove_node(v)
    return cycles
      

#### crisp with bounded forks ####

