from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import networkx as nx

from . import probs as PROBS


SOURCE_NODE = "__SOURCE__"


def edge_key(u: Any, v: Any) -> tuple:
    return tuple(sorted((u, v), key=repr))


@dataclass(frozen=True)
class ChainSummary:
    endpoints: tuple[Any, Any]
    nodes: tuple[Any, ...]
    edges: tuple[tuple[Any, Any], ...]
    internal_nodes: tuple[Any, ...]
    length: float
    weight: float
    edge_weight: float
    original_nodes: frozenset[Any] = field(default_factory=frozenset)
    original_edges: frozenset[tuple[Any, Any]] = field(default_factory=frozenset)

    @property
    def total_weight(self) -> float:
        return self.weight + self.edge_weight


@dataclass
class EdgeReliabilityDecomposition:
    source_graph: nx.Graph
    source: Any
    bridges: tuple[tuple[Any, Any], ...]
    two_edge_components: tuple[frozenset[Any], ...]
    bridge_tree: nx.Graph
    structure_graph: nx.Graph
    three_edge_macro_graph: nx.Graph
    regular_chains: tuple[ChainSummary, ...]
    generalized_chains: tuple[ChainSummary, ...]
    component_sources: dict[int, tuple[Any, ...]]

    def td_saidi(self, *, graph: nx.Graph | None = None, source: Any | None = None, max_fail: int = 2):
        """Compute a truncated TD SAIDI polynomial on one decomposition graph.

        This helper intentionally supports node weights only. If a graph contains
        nonzero `edge_weight`, exact reliability should be handled by a dedicated
        weighted contraction, not silently by `GraphRel`.
        """
        from .graph_rel import GraphRel

        target_graph = self.source_graph if graph is None else graph
        target_source = self.source if source is None else source
        if isinstance(target_graph, nx.MultiGraph):
            raise TypeError("decomposition TD helper does not support MultiGraph")
        if any(float(data.get("edge_weight", 0.0)) != 0.0 for _, _, data in target_graph.edges(data=True)):
            raise NotImplementedError("TD helper currently supports node weights only, not edge_weight")
        edge_probs = {edge: PROBS.Poly([0, 1]) for edge in target_graph.edges}
        node_weights = {node: float(data.get("weight", 0.0)) for node, data in target_graph.nodes(data=True)}
        graph_rel = GraphRel(
            target_graph,
            nodes_weight=node_weights,
            edges_prob=edge_probs,
            max_fail=max_fail,
            sources=[target_source],
        )
        return graph_rel.calc_rel(["saidi"])["saidi"]

    def tree_td_saidi(self, *, max_fail: int = 2):
        source_component = self.source_graph.nodes[self.source]["component_id"]
        return self.td_saidi(graph=self.bridge_tree, source=source_component, max_fail=max_fail)


def decompose_graph_rel(
    graph_rel,
    *,
    include_generalized_chains: bool = True,
) -> EdgeReliabilityDecomposition:
    if isinstance(graph_rel.original_graph, nx.MultiGraph):
        raise TypeError("GraphRel decomposition currently supports nx.Graph input only, not MultiGraph")

    source_graph = _simple_source_graph(graph_rel)
    source = graph_rel.source
    bridges = tuple(sorted((edge_key(u, v) for u, v in nx.bridges(source_graph)), key=repr))
    component_sets, node_to_component = _two_edge_components(source_graph, bridges)
    _set_component_metadata(source_graph, component_sets, node_to_component)
    bridge_tree = _build_bridge_tree(source_graph, bridges, component_sets, node_to_component)
    component_sources = _component_sources(source_graph, bridges, node_to_component, source)
    structure_graph, regular_chains = _build_structure_graph(
        source_graph,
        bridges,
        component_sets,
        component_sources,
    )
    if include_generalized_chains:
        three_edge_macro_graph = _contract_three_edge_components(structure_graph, sources=[source])
        generalized_chains = extract_chains(three_edge_macro_graph, sources=[source])
    else:
        three_edge_macro_graph = nx.Graph()
        generalized_chains = []

    return EdgeReliabilityDecomposition(
        source_graph=source_graph,
        source=source,
        bridges=bridges,
        two_edge_components=tuple(frozenset(comp) for comp in component_sets),
        bridge_tree=bridge_tree,
        structure_graph=structure_graph,
        three_edge_macro_graph=three_edge_macro_graph,
        regular_chains=tuple(regular_chains),
        generalized_chains=tuple(generalized_chains),
        component_sources=component_sources,
    )


def _simple_source_graph(graph_rel) -> nx.Graph:
    graph = graph_rel.graph.copy()
    graph.remove_edges_from(nx.selfloop_edges(graph))

    edge_counts: dict[tuple, int] = {}
    for u, v, _ in graph.edges(keys=True):
        key = edge_key(u, v)
        edge_counts[key] = edge_counts.get(key, 0) + 1
    parallel_edges = [edge for edge, count in edge_counts.items() if count > 1]
    if parallel_edges:
        raise TypeError(
            "GraphRel decomposition does not support parallel edges after source contraction; "
            f"examples: {parallel_edges[:5]}"
        )

    simple = nx.Graph()
    for node, data in graph.nodes(data=True):
        node_data = data.copy()
        node_data.setdefault("original_nodes", frozenset(data.get("component", {node})))
        node_data["weight"] = float(node_data.get("weight", 0.0))
        simple.add_node(node, **node_data)
    for u, v, data in graph.edges(data=True):
        edge_data = data.copy()
        edge_data.setdefault("original_edges", frozenset({edge_key(u, v)}))
        edge_data.setdefault("length", edge_data.get("length_m", 1.0))
        edge_data.setdefault("edge_weight", 0.0)
        simple.add_edge(u, v, **edge_data)
    return simple


def _two_edge_components(
    graph: nx.Graph,
    bridges: Iterable[tuple[Any, Any]],
) -> tuple[list[set[Any]], dict[Any, int]]:
    without_bridges = graph.copy()
    without_bridges.remove_edges_from(list(bridges))
    components = [set(comp) for comp in nx.connected_components(without_bridges)]
    node_to_component = {}
    for idx, comp in enumerate(components):
        for node in comp:
            node_to_component[node] = idx
    return components, node_to_component


def _set_component_metadata(
    graph: nx.Graph,
    components: list[set[Any]],
    node_to_component: dict[Any, int],
) -> None:
    for node in graph.nodes:
        graph.nodes[node]["component_id"] = node_to_component[node]
    for idx, comp in enumerate(components):
        weight = sum(float(graph.nodes[node].get("weight", 0.0)) for node in comp)
        for node in comp:
            graph.nodes[node]["component_weight"] = weight


def _build_bridge_tree(
    graph: nx.Graph,
    bridges: Iterable[tuple[Any, Any]],
    components: list[set[Any]],
    node_to_component: dict[Any, int],
) -> nx.Graph:
    bridge_tree = nx.Graph()
    for idx, comp in enumerate(components):
        original_nodes = _collect_original_nodes(graph, comp)
        bridge_tree.add_node(
            idx,
            weight=sum(float(graph.nodes[node].get("weight", 0.0)) for node in comp),
            original_nodes=frozenset(original_nodes),
            members=frozenset(comp),
        )
    for u, v in bridges:
        cu = node_to_component[u]
        cv = node_to_component[v]
        if cu == cv:
            continue
        data = graph.edges[u, v].copy()
        data["original_edge"] = edge_key(u, v)
        data["original_edges"] = frozenset({edge_key(u, v)})
        bridge_tree.add_edge(cu, cv, **data)
    return bridge_tree


def _component_sources(
    graph: nx.Graph,
    bridges: Iterable[tuple[Any, Any]],
    node_to_component: dict[Any, int],
    source: Any,
) -> dict[int, tuple[Any, ...]]:
    sources: dict[int, set[Any]] = defaultdict_set()
    sources[node_to_component[source]].add(source)
    for u, v in bridges:
        sources[node_to_component[u]].add(u)
        sources[node_to_component[v]].add(v)
    return {component: tuple(sorted(nodes, key=repr)) for component, nodes in sources.items()}


def defaultdict_set():
    from collections import defaultdict

    return defaultdict(set)


def _build_structure_graph(
    graph: nx.Graph,
    bridges: Iterable[tuple[Any, Any]],
    components: list[set[Any]],
    component_sources: dict[int, tuple[Any, ...]],
) -> tuple[nx.Graph, tuple[ChainSummary, ...]]:
    structure = nx.Graph()
    all_chains: list[ChainSummary] = []
    for idx, comp in enumerate(components):
        subgraph = graph.subgraph(comp).copy()
        if subgraph.number_of_edges() == 0:
            continue
        sources = component_sources.get(idx, ())
        chains = extract_chains(subgraph, sources=sources)
        all_chains.extend(chains)
        for chain in chains:
            u, v = chain.endpoints
            for node in (u, v):
                if node not in structure:
                    data = graph.nodes[node].copy()
                    data["weight"] = float(data.get("weight", 0.0))
                    data.setdefault("original_nodes", frozenset({node}))
                    structure.add_node(node, **data)
            if u == v:
                continue
            if structure.has_edge(u, v):
                existing = structure.edges[u, v]
                existing["parallel_chain_count"] = int(existing.get("parallel_chain_count", 1)) + 1
                existing["length"] = float(existing.get("length", 0.0)) + chain.length
                existing["edge_weight"] = float(existing.get("edge_weight", 0.0)) + chain.total_weight
                existing["original_nodes"] = frozenset(set(existing.get("original_nodes", set())) | set(chain.original_nodes))
                existing["original_edges"] = frozenset(set(existing.get("original_edges", set())) | set(chain.original_edges))
                existing["parallel_chain_nodes"] = tuple(existing.get("parallel_chain_nodes", ())) + (chain.nodes,)
                existing["parallel_chain_edges"] = tuple(existing.get("parallel_chain_edges", ())) + (chain.edges,)
                continue
            structure.add_edge(
                u,
                v,
                length=chain.length,
                edge_weight=chain.total_weight,
                parallel_chain_count=1,
                chain_nodes=chain.nodes,
                chain_edges=chain.edges,
                original_nodes=chain.original_nodes,
                original_edges=chain.original_edges,
                is_bridge=False,
            )
    return structure, tuple(all_chains)


def _contract_three_edge_components(graph: nx.Graph, sources: Iterable[Any]) -> nx.Graph:
    if graph.number_of_nodes() == 0:
        return nx.Graph()

    source_set = set(sources)
    components = [set(comp) for comp in nx.k_edge_components(graph, k=3)]
    node_to_component = {}
    for idx, comp in enumerate(components):
        for node in comp:
            node_to_component[node] = idx

    macro = nx.Graph()
    for idx, comp in enumerate(components):
        original_nodes = _collect_original_nodes(graph, comp)
        internal_edge_weight = sum(
            float(data.get("edge_weight", 0.0))
            for u, v, data in graph.edges(comp, data=True)
            if u in comp and v in comp
        )
        macro_node = _component_node_name("C3", idx, comp, source_set)
        macro.add_node(
            macro_node,
            weight=sum(float(graph.nodes[node].get("weight", 0.0)) for node in comp) + internal_edge_weight,
            original_nodes=frozenset(original_nodes),
            members=frozenset(comp),
            contains_source=bool(comp & source_set),
        )
        for node in comp:
            node_to_component[node] = macro_node

    for u, v, data in graph.edges(data=True):
        cu = node_to_component[u]
        cv = node_to_component[v]
        if cu == cv:
            continue
        if macro.has_edge(cu, cv):
            raise TypeError(
                "3-edge macro graph would contain parallel edges; "
                "MultiGraph decomposition is intentionally unsupported"
            )
        edge_data = data.copy()
        edge_data.setdefault("edge_weight", 0.0)
        edge_data.setdefault("length", data.get("length_m", 1.0))
        macro.add_edge(cu, cv, **edge_data)
    return macro


def _component_node_name(prefix: str, idx: int, comp: set[Any], sources: set[Any]) -> Any:
    source_hits = comp & sources
    if len(comp) == 1:
        return next(iter(comp))
    if source_hits:
        return next(iter(source_hits))
    return f"{prefix}_{idx}"


def extract_chains(
    graph: nx.Graph,
    *,
    sources: Iterable[Any] = (),
    node_weight_attr: str = "weight",
    edge_weight_attr: str = "edge_weight",
    length_attr: str = "length",
) -> list[ChainSummary]:
    if isinstance(graph, nx.MultiGraph):
        raise TypeError("chain extraction does not support MultiGraph")
    if graph.number_of_nodes() == 0:
        return []

    source_set = set(sources)
    boundary_nodes = {node for node, degree in graph.degree if degree != 2} | (source_set & set(graph.nodes))
    for comp in nx.connected_components(graph):
        if boundary_nodes.isdisjoint(comp):
            boundary_nodes.add(next(iter(comp)))

    visited_edges: set[tuple[Any, Any]] = set()
    chains: list[ChainSummary] = []
    for start in sorted(boundary_nodes, key=repr):
        for neighbor in sorted(graph.neighbors(start), key=repr):
            key = edge_key(start, neighbor)
            if key in visited_edges:
                continue
            chains.append(
                _walk_chain(
                    graph,
                    start=start,
                    next_node=neighbor,
                    boundary_nodes=boundary_nodes,
                    visited_edges=visited_edges,
                    node_weight_attr=node_weight_attr,
                    edge_weight_attr=edge_weight_attr,
                    length_attr=length_attr,
                )
            )
    return chains


def _walk_chain(
    graph: nx.Graph,
    *,
    start: Any,
    next_node: Any,
    boundary_nodes: set[Any],
    visited_edges: set[tuple[Any, Any]],
    node_weight_attr: str,
    edge_weight_attr: str,
    length_attr: str,
) -> ChainSummary:
    nodes = [start, next_node]
    edges = [edge_key(start, next_node)]
    visited_edges.add(edge_key(start, next_node))
    prev = start
    current = next_node

    while current not in boundary_nodes:
        candidates = [
            neighbor
            for neighbor in graph.neighbors(current)
            if neighbor != prev and edge_key(current, neighbor) not in visited_edges
        ]
        if not candidates:
            break
        nxt = sorted(candidates, key=repr)[0]
        visited_edges.add(edge_key(current, nxt))
        edges.append(edge_key(current, nxt))
        nodes.append(nxt)
        prev, current = current, nxt

    if current == start and len(nodes) > 1:
        internal_nodes = tuple(nodes[1:-1])
    else:
        internal_nodes = tuple(nodes[1:-1])
    length = sum(float(graph.edges[edge].get(length_attr, graph.edges[edge].get("length_m", 1.0))) for edge in edges)
    edge_weight = sum(float(graph.edges[edge].get(edge_weight_attr, 0.0)) for edge in edges)
    node_weight = sum(float(graph.nodes[node].get(node_weight_attr, 0.0)) for node in internal_nodes)
    original_nodes = set()
    for node in nodes:
        original_nodes.update(graph.nodes[node].get("original_nodes", {node}))
    original_edges = set()
    for edge in edges:
        original_edges.update(graph.edges[edge].get("original_edges", {edge}))
    return ChainSummary(
        endpoints=(start, current),
        nodes=tuple(nodes),
        edges=tuple(edges),
        internal_nodes=internal_nodes,
        length=length,
        weight=node_weight,
        edge_weight=edge_weight,
        original_nodes=frozenset(original_nodes),
        original_edges=frozenset(original_edges),
    )


def _collect_original_nodes(graph: nx.Graph, nodes: Iterable[Any]) -> set[Any]:
    original_nodes = set()
    for node in nodes:
        original_nodes.update(graph.nodes[node].get("original_nodes", {node}))
    return original_nodes
