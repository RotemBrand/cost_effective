from __future__ import annotations

from typing import Iterable

import networkx as nx
import numpy as np


def add_synthetic_switches(
    graph: nx.Graph,
    *,
    n_switches: int,
    seed: int | None = None,
    length_attr: str = "length_m",
    pos_attr: str = "pos",
) -> nx.Graph:
    """Return a copy of ``graph`` with synthetic switch metadata.

    The topology is not changed. Non-tree edges are marked as normally open tie
    switches. Then a random subset of normally closed edges is marked as
    sectionalizing switches so the total number of ``is_switch`` edges is close
    to ``n_switches``. If ``n_switches`` is smaller than the number of ties, all
    ties remain switches.
    """
    if n_switches < 0:
        raise ValueError("n_switches must be non-negative")
    if isinstance(graph, nx.MultiGraph):
        raise TypeError("add_synthetic_switches currently supports nx.Graph, not nx.MultiGraph")

    graph_copy = nx.Graph(graph).copy()
    ensure_edge_lengths(graph_copy, length_attr=length_attr, pos_attr=pos_attr)

    tree_edge_keys = {
        _edge_key(edge)
        for edge in nx.minimum_spanning_tree(graph_copy, weight=length_attr).edges
    }
    tie_edges = [
        (u, v)
        for u, v in graph_copy.edges
        if _edge_key((u, v)) not in tree_edge_keys
    ]
    tie_edge_keys = {_edge_key(edge) for edge in tie_edges}

    for u, v, data in graph_copy.edges(data=True):
        is_tie = _edge_key((u, v)) in tie_edge_keys
        data["is_tie"] = is_tie
        data["is_switch"] = is_tie
        data["normally_closed"] = not is_tie
        data["synthetic_switch"] = False

    closed_edges = [
        (u, v)
        for u, v in graph_copy.edges
        if _edge_key((u, v)) not in tie_edge_keys
    ]
    n_closed_switches = min(max(n_switches - len(tie_edges), 0), len(closed_edges))
    rng = np.random.default_rng(seed)
    if n_closed_switches:
        selected_idx = rng.choice(len(closed_edges), size=n_closed_switches, replace=False)
        for idx in selected_idx:
            u, v = closed_edges[int(idx)]
            graph_copy.edges[u, v]["is_switch"] = True
            graph_copy.edges[u, v]["synthetic_switch"] = True

    graph_copy.graph.update(graph.graph)
    graph_copy.graph["synthetic_switches"] = True
    graph_copy.graph["requested_n_switches"] = int(n_switches)
    graph_copy.graph["n_tie_switches"] = len(tie_edges)
    graph_copy.graph["n_closed_switches"] = int(n_closed_switches)
    graph_copy.graph["n_switches"] = sum(
        1 for _, _, data in graph_copy.edges(data=True) if data.get("is_switch", False)
    )
    return graph_copy


def add_synthetic_switches_like(
    target_graph: nx.Graph,
    reference_graph: nx.Graph,
    *,
    seed: int | None = None,
    length_attr: str = "length_m",
    pos_attr: str = "pos",
) -> nx.Graph:
    """Add synthetic switches using the total switch count of ``reference_graph``."""
    return add_synthetic_switches(
        target_graph,
        n_switches=count_switch_edges(reference_graph),
        seed=seed,
        length_attr=length_attr,
        pos_attr=pos_attr,
    )


def count_switch_edges(graph: nx.Graph) -> int:
    return sum(1 for _, _, data in graph.edges(data=True) if data.get("is_switch", False))


def count_tie_edges(graph: nx.Graph) -> int:
    return sum(1 for _, _, data in graph.edges(data=True) if data.get("is_tie", False))


def ensure_edge_lengths(
    graph: nx.Graph,
    *,
    length_attr: str = "length_m",
    pos_attr: str = "pos",
) -> None:
    """Set edge ``length_attr`` in-place from existing attrs or node positions."""
    for u, v, data in graph.edges(data=True):
        length = _edge_length(graph, u, v, data, length_attr=length_attr, pos_attr=pos_attr)
        data[length_attr] = length
        data.setdefault("length", length)


def _edge_key(edge: Iterable) -> tuple:
    return tuple(sorted(tuple(edge)[:2]))


def _edge_length(
    graph: nx.Graph,
    u,
    v,
    data: dict,
    *,
    length_attr: str,
    pos_attr: str,
) -> float:
    for attr in (length_attr, "length_m", "length", "weight"):
        value = data.get(attr)
        if value is not None:
            try:
                return float(value)
            except (TypeError, ValueError):
                pass
    pos = nx.get_node_attributes(graph, pos_attr)
    if u not in pos or v not in pos:
        return 1.0
    return float(np.linalg.norm(np.asarray(pos[u], dtype=float) - np.asarray(pos[v], dtype=float)))
