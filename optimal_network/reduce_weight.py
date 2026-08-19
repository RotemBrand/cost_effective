import numpy as np
import networkx as nx
from indexes import get_skeleton_graph
from itertools import product

def optimize_rel_weight_ratio(optimal_graph: nx.Graph, max_risk_gain: float=0.5, max_changes: int=300, source=None, debug: bool=False) -> nx.Graph:    
    new_graph = optimal_graph.copy()
    weight = _total_weight(optimal_graph)
    risk = _total_inter_risk(optimal_graph)
    risk_gain = 0

    for it in range(max_changes):
        # get optimal change
        changes_data = _get_all_changes(new_graph, risk, max_risk_gain - risk_gain, source=source)
        if len(changes_data) == 0:
            break
        changes_data = sorted(changes_data, key=lambda row: -row['score'])

        optimal_row = max(changes_data, key=lambda row: row['score'])


        # make change
        removed_attrs = dict(new_graph.edges[optimal_row['fork_neigh'], optimal_row['fork']])
        new_graph.add_edge(optimal_row['fork_neigh'], optimal_row['new_fork'])
        new_graph.remove_edge(optimal_row['fork_neigh'], optimal_row['fork'])
        
        risk_gain = (_total_inter_risk(new_graph) - risk) / risk
        weight_gain = (_total_weight(new_graph) - weight) / weight

        # debug
        if debug:
            print(f"{it}: {weight_gain=:.3f}, {risk_gain=:.3f}")
        if risk_gain > max_risk_gain:
            new_graph.remove_edge(optimal_row['fork_neigh'], optimal_row['new_fork'])
            new_graph.add_edge(optimal_row['fork_neigh'], optimal_row['fork'], **removed_attrs)
            break

    new_weight = _total_weight(new_graph)
    new_risk = _total_inter_risk(new_graph)
    print(f"weight_gain = {(new_weight - weight) / weight:.2f}, risk_gain = {(new_risk - risk) / risk}")
    return new_graph

    


def _get_all_changes(optimal_graph: nx.Graph, original_risk: float, max_risk_gain: float, source=None) -> list:
    def get_weight(node1, node2) -> float:
        pos1 = np.asarray(optimal_graph.nodes[node1]["pos"], dtype=float)
        pos2 = np.asarray(optimal_graph.nodes[node2]["pos"], dtype=float)
        return np.linalg.norm(pos1 - pos2)
    # optimal graph stats
    total_weight = _total_weight(optimal_graph)
    inter_risk = _total_inter_risk(optimal_graph)
    
    # get fork chains data
    fork_chains = _forks_chains_data(optimal_graph, source=source)
    changes_data = [] # fork, fork_neigh, new_fork, score
    for fork, f_data in fork_chains.items():
        for c_data1, c_data2 in product(f_data, f_data):
            chain1, chain2 = c_data1["chain"], c_data2["chain"]
            if chain1 == chain2 or c_data1["len"] == 0 or c_data2["len"] == 0:
                continue
            fork_neigh = c_data1["nodes_order"][1]
            chain1_len = c_data1["len"]
            chain2_len = c_data2["len"]
            for new_fork_idx, new_fork in enumerate(c_data2["nodes_order"]):
                if new_fork in {fork, fork_neigh} or optimal_graph.has_edge(fork_neigh, new_fork):
                    continue
                weight_gain = (get_weight(fork_neigh, new_fork) - get_weight(fork_neigh, fork)) / total_weight
                if weight_gain >= 0:
                    continue
               
                risk_gain = (
                    _inter_risk(chain2_len - new_fork_idx) - _inter_risk(chain2_len) + 
                    _inter_risk(chain1_len + new_fork_idx) - _inter_risk(chain1_len)
                ) / original_risk
                if risk_gain > max_risk_gain:
                    continue

                score = -weight_gain / risk_gain if risk_gain != 0 else np.inf
                changes_data.append({
                    'fork': fork,
                    'fork_neigh': fork_neigh,
                    'new_fork': new_fork,
                    'weight_gain': weight_gain,
                    'risk_gain': risk_gain,
                    'score': score
                })
    return changes_data



def _forks_chains_data(optimal_graph: nx.Graph, source=None) -> dict:
    sources = [source] if source is not None else None
    strc = get_skeleton_graph(optimal_graph, sources=sources)
    forks_chains = {}
    for fork in strc.nodes:
        if (source is not None) and fork == source:
            continue
        forks_chains[fork] = []
        for c_and_c_data in strc.edges(fork, data=True, keys=True):
            chain = _key_chain(c_and_c_data[:3])
            c_data = c_and_c_data[3]
            nodes_order = list(nx.dfs_preorder_nodes(c_data["subgraph"], fork))
            forks_chains[fork].append({'chain': chain, "nodes_order": nodes_order, "len": c_data["length"]})
    return forks_chains

def _key_chain(c: tuple) -> tuple:
    return (*sorted(c[:2]), c[2])

def _inter_risk(l: int) -> float:
    n = l - 1
    return 1 / 6 * (n ** 3 + 3 * n**2 + 2 * n)

def _total_weight(G: nx.Graph) -> float:
    pos = nx.get_node_attributes(G, "pos")
    return sum([
        np.linalg.norm(np.asarray(pos[v], dtype=float) - np.asarray(pos[u], dtype=float))
        for v, u in G.edges
    ])

def _total_inter_risk(G: nx.Graph) -> float:
    strc = get_skeleton_graph(G)
    return sum(map(_inter_risk, nx.get_edge_attributes(strc, "length").values()))

