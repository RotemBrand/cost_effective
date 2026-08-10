from __future__ import annotations

from functools import reduce
from operator import mul
from typing import Iterable

import networkx as nx
import numpy as np


def key_edge(edge: tuple) -> tuple:
    return tuple(sorted(edge[:2]))


def contract_switch_sections(
    graph: nx.Graph,
    sources: Iterable,
    *,
    node_weight_attr: str = "weight",
    edge_weight_attr: str = "edge_weight",
    tie_attr: str = "is_tie",
    switch_attr: str | None = "is_switch",
    prob_attr: str = "prob",
    length_attr: str = "length",
    pos_attr: str = "pos",
) -> nx.Graph:
    """Contract non-tie degree-2 chains into switchable section edges.

    Sources, non-degree-2 nodes, and endpoints of switch or tie edges remain explicit.
    Internal node weights are moved to the contracted section edge's
    ``edge_weight``. Tie edges are copied explicitly and are never swallowed into
    an ordinary section.
    """
    if isinstance(graph, nx.MultiGraph):
        raise TypeError("contract_switch_sections currently supports nx.Graph, not nx.MultiGraph")

    sources = set(sources)
    if not sources <= set(graph.nodes):
        missing = sources - set(graph.nodes)
        raise ValueError(f"sources are not in graph: {missing}")

    graph_copy = nx.Graph(graph)
    _set_default_node_weights(graph_copy, node_weight_attr)
    _set_default_edge_attrs(graph_copy, edge_weight_attr, tie_attr)

    boundaries = _section_boundaries(
        graph_copy,
        sources=sources,
        tie_attr=tie_attr,
        switch_attr=switch_attr,
    )
    if len(boundaries) == 0 and graph_copy.number_of_nodes() > 0:
        boundaries.add(next(iter(graph_copy.nodes)))

    section_graph = nx.MultiGraph()
    for node in boundaries:
        data = graph_copy.nodes[node].copy()
        section_graph.add_node(node, **data)

    visited_non_tie_edges: set[tuple] = set()

    for u, v, data in graph_copy.edges(data=True):
        if data.get(tie_attr, False):
            section_graph.add_edge(
                u,
                v,
                **_tie_edge_attrs(graph_copy, (u, v), edge_weight_attr, tie_attr),
            )

    for boundary in list(boundaries):
        for neighbor in graph_copy.neighbors(boundary):
            edge = key_edge((boundary, neighbor))
            if edge in visited_non_tie_edges:
                continue
            if graph_copy.edges[boundary, neighbor].get(tie_attr, False):
                continue
            path_nodes, path_edges = _walk_section_chain(
                graph_copy,
                start=boundary,
                next_node=neighbor,
                boundaries=boundaries,
                tie_attr=tie_attr,
            )
            visited_non_tie_edges.update(path_edges)
            _add_section_edge(
                section_graph,
                graph_copy,
                path_nodes,
                path_edges,
                node_weight_attr=node_weight_attr,
                edge_weight_attr=edge_weight_attr,
                tie_attr=tie_attr,
                prob_attr=prob_attr,
                length_attr=length_attr,
                pos_attr=pos_attr,
            )

    for u, v, data in graph_copy.edges(data=True):
        edge = key_edge((u, v))
        if data.get(tie_attr, False) or edge in visited_non_tie_edges:
            continue
        section_graph.add_node(u, **graph_copy.nodes[u].copy())
        section_graph.add_node(v, **graph_copy.nodes[v].copy())
        section_graph.add_edge(
            u,
            v,
            **_section_edge_attrs(
                graph_copy,
                [u, v],
                [edge],
                node_weight_attr=node_weight_attr,
                edge_weight_attr=edge_weight_attr,
                tie_attr=tie_attr,
                prob_attr=prob_attr,
                length_attr=length_attr,
                pos_attr=pos_attr,
            ),
        )
        visited_non_tie_edges.add(edge)

    section_graph.graph.update(graph.graph)
    section_graph.graph["sources"] = list(sources)
    section_graph.graph["section_contracted"] = True
    return section_graph


def _set_default_node_weights(graph: nx.Graph, node_weight_attr: str) -> None:
    for node in graph.nodes:
        graph.nodes[node].setdefault(node_weight_attr, 1.0)


def _set_default_edge_attrs(graph: nx.Graph, edge_weight_attr: str, tie_attr: str) -> None:
    for _, _, data in graph.edges(data=True):
        data.setdefault(edge_weight_attr, 0.0)
        data.setdefault(tie_attr, False)


def _section_boundaries(
    graph: nx.Graph,
    sources: set,
    tie_attr: str,
    switch_attr: str | None,
) -> set:
    boundaries = set(sources)
    boundaries.update(node for node, degree in graph.degree if degree != 2)
    for u, v, data in graph.edges(data=True):
        if data.get(tie_attr, False) or (switch_attr is not None and data.get(switch_attr, False)):
            boundaries.add(u)
            boundaries.add(v)
    return boundaries


def _walk_section_chain(
    graph: nx.Graph,
    *,
    start,
    next_node,
    boundaries: set,
    tie_attr: str,
) -> tuple[list, list[tuple]]:
    path_nodes = [start, next_node]
    path_edges = [key_edge((start, next_node))]
    prev = start
    current = next_node

    while current not in boundaries:
        candidates = [
            n
            for n in graph.neighbors(current)
            if n != prev and not graph.edges[current, n].get(tie_attr, False)
        ]
        if len(candidates) == 0:
            break
        if len(candidates) > 1:
            break
        nxt = candidates[0]
        path_nodes.append(nxt)
        path_edges.append(key_edge((current, nxt)))
        prev, current = current, nxt

    return path_nodes, path_edges


def _add_section_edge(
    section_graph: nx.Graph,
    graph: nx.Graph,
    path_nodes: list,
    path_edges: list[tuple],
    *,
    node_weight_attr: str,
    edge_weight_attr: str,
    tie_attr: str,
    prob_attr: str,
    length_attr: str,
    pos_attr: str,
) -> None:
    end_a = path_nodes[0]
    end_b = path_nodes[-1]
    for node in (end_a, end_b):
        if node not in section_graph:
            section_graph.add_node(node, **graph.nodes[node].copy())
    section_graph.add_edge(
        end_a,
        end_b,
        **_section_edge_attrs(
            graph,
            path_nodes,
            path_edges,
            node_weight_attr=node_weight_attr,
            edge_weight_attr=edge_weight_attr,
            tie_attr=tie_attr,
            prob_attr=prob_attr,
            length_attr=length_attr,
            pos_attr=pos_attr,
        ),
    )


def _section_edge_attrs(
    graph: nx.Graph,
    path_nodes: list,
    path_edges: list[tuple],
    *,
    node_weight_attr: str,
    edge_weight_attr: str,
    tie_attr: str,
    prob_attr: str,
    length_attr: str,
    pos_attr: str,
) -> dict:
    internal_nodes = path_nodes[1:-1]
    edge_attrs = [graph.edges[e] for e in path_edges]
    attrs = {
        "is_section": True,
        tie_attr: False,
        "original_nodes": list(path_nodes),
        "internal_nodes": list(internal_nodes),
        "original_edges": list(path_edges),
        edge_weight_attr: (
            sum(float(graph.nodes[n].get(node_weight_attr, 0.0)) for n in internal_nodes)
            + sum(float(data.get(edge_weight_attr, 0.0)) for data in edge_attrs)
        ),
        length_attr: _path_length(graph, path_edges, length_attr=length_attr, pos_attr=pos_attr),
    }
    if all(prob_attr in data for data in edge_attrs):
        attrs[prob_attr] = _series_failure_probability([data[prob_attr] for data in edge_attrs])
    return attrs


def _tie_edge_attrs(
    graph: nx.Graph,
    edge: tuple,
    edge_weight_attr: str,
    tie_attr: str,
) -> dict:
    data = graph.edges[edge].copy()
    data.setdefault(edge_weight_attr, 0.0)
    data[tie_attr] = True
    data["is_section"] = False
    data["original_nodes"] = list(edge)
    data["original_edges"] = [key_edge(edge)]
    return data


def _path_length(
    graph: nx.Graph,
    path_edges: list[tuple],
    *,
    length_attr: str,
    pos_attr: str,
) -> float:
    return sum(_edge_length(graph, edge, length_attr=length_attr, pos_attr=pos_attr) for edge in path_edges)


def _edge_length(graph: nx.Graph, edge: tuple, *, length_attr: str, pos_attr: str) -> float:
    data = graph.edges[edge]
    if length_attr in data and data[length_attr] is not None:
        return float(data[length_attr])
    if "length_m" in data and data["length_m"] is not None:
        return float(data["length_m"])
    pos = nx.get_node_attributes(graph, pos_attr)
    if edge[0] in pos and edge[1] in pos:
        return float(np.linalg.norm(np.asarray(pos[edge[0]]) - np.asarray(pos[edge[1]])))
    return 1.0


def _series_failure_probability(probabilities: list):
    success_probability = reduce(mul, [1 - p for p in probabilities], 1)
    return 1 - success_probability
