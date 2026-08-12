from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys
from typing import Any, Iterable, Literal

import networkx as nx

from . import probs as PROBS


SOURCE_NODE = "__SOURCE__"
ThreeEdgeComponentMethod = Literal["networkx", "projection"]


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
    generalized_component_method: str = "networkx"

    @property
    def total_weight(self) -> float:
        return _total_graph_weight(self.source_graph)

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

    def switch_risk_terms(
        self,
        *,
        edge_probs: dict | None = None,
        mean_edge_failure_prob: float | None = None,
        length_attr: str = "length",
        output: Literal["float", "poly"] = "float",
    ) -> "SwitchRiskTerms":
        """Return tree and chain leading risk terms for a switch-section graph.

        `output="float"` evaluates the terms using finite edge probabilities.
        `output="poly"` returns leading polynomials in a mean-probability
        parameter, with coefficients through p^2.
        """
        if output not in ("float", "poly"):
            raise ValueError("output must be 'float' or 'poly'")
        prepared = _prepare_edge_probabilities(
            self.source_graph,
            edge_probs=edge_probs,
            mean_edge_failure_prob=mean_edge_failure_prob,
            length_attr=length_attr,
            output=output,
        )
        tree = _tree_risk(
            self.bridge_tree,
            self.source_graph.nodes[self.source]["component_id"],
            prepared,
            total_weight=self.total_weight,
            output=output,
        )
        section = _nonbridge_section_edge_weight_risk(
            self.source_graph,
            self.bridges,
            prepared,
            total_weight=self.total_weight,
            output=output,
        )
        internal = _chain_two_cut_risk(
            self.regular_chains,
            self.source_graph,
            prepared,
            total_weight=self.total_weight,
            output=output,
        )
        structural = _chain_two_cut_risk(
            self.generalized_chains,
            self.three_edge_macro_graph,
            prepared,
            total_weight=self.total_weight,
            output=output,
        )
        return SwitchRiskTerms(
            total_weight=self.total_weight,
            tree=tree,
            nonbridge_section=section,
            internal=internal,
            structural=structural,
            output=output,
            mean_edge_failure_prob=mean_edge_failure_prob,
            edge_failure_rate=prepared.edge_failure_rate,
            mean_actual_edge_probability=prepared.mean_actual_edge_probability,
        )


@dataclass(frozen=True)
class SwitchRiskTerms:
    total_weight: float
    tree: Any
    nonbridge_section: Any
    internal: Any
    structural: Any
    output: str
    mean_edge_failure_prob: float | None = None
    edge_failure_rate: float | None = None
    mean_actual_edge_probability: float | None = None

    @property
    def total(self):
        return self.tree + self.nonbridge_section + self.internal + self.structural


@dataclass(frozen=True)
class _PreparedEdgeProbabilities:
    q: dict[tuple[Any, Any], float]
    alpha: dict[tuple[Any, Any], float]
    edge_failure_rate: float | None
    mean_actual_edge_probability: float | None


def decompose_graph_rel(
    graph_rel,
    *,
    include_generalized_chains: bool = True,
    generalized_component_method: ThreeEdgeComponentMethod = "networkx",
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
        three_edge_macro_graph = _contract_three_edge_components(
            structure_graph,
            sources=[source],
            chains=regular_chains,
            component_method=generalized_component_method,
        )
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
        generalized_component_method=generalized_component_method,
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
        internal_edge_weight = sum(
            float(data.get("edge_weight", 0.0))
            for u, v, data in graph.edges(comp, data=True)
            if u in comp and v in comp
        )
        bridge_tree.add_node(
            idx,
            weight=sum(float(graph.nodes[node].get("weight", 0.0)) for node in comp) + internal_edge_weight,
            node_weight=sum(float(graph.nodes[node].get("weight", 0.0)) for node in comp),
            internal_edge_weight=internal_edge_weight,
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


def _total_graph_weight(graph: nx.Graph) -> float:
    return sum(float(data.get("weight", 0.0)) for _, data in graph.nodes(data=True)) + sum(
        float(data.get("edge_weight", 0.0)) for _, _, data in graph.edges(data=True)
    )


def _prepare_edge_probabilities(
    graph: nx.Graph,
    *,
    edge_probs: dict | None,
    mean_edge_failure_prob: float | None,
    length_attr: str,
    output: Literal["float", "poly"],
) -> _PreparedEdgeProbabilities:
    if edge_probs is not None and mean_edge_failure_prob is not None:
        raise ValueError("pass either edge_probs or mean_edge_failure_prob, not both")

    edge_failure_rate = None
    mean_actual = None
    if mean_edge_failure_prob is not None:
        from .utilities import edge_probs_by_length

        q_raw, edge_failure_rate = edge_probs_by_length(
            graph,
            p=mean_edge_failure_prob,
            mode="mean",
            length_attr=length_attr,
        )
        q = {_canonical_edge_key(edge): float(prob) for edge, prob in q_raw.items()}
        if mean_edge_failure_prob == 0:
            alpha = {edge: 0.0 for edge in q}
        else:
            alpha = {edge: prob / float(mean_edge_failure_prob) for edge, prob in q.items()}
        mean_actual = sum(q.values()) / len(q) if q else 0.0
        return _PreparedEdgeProbabilities(
            q=q,
            alpha=alpha,
            edge_failure_rate=edge_failure_rate,
            mean_actual_edge_probability=mean_actual,
        )

    if edge_probs is None:
        edge_probs = {
            edge_key(u, v): data.get("prob", PROBS.Poly([0, 1]))
            for u, v, data in graph.edges(data=True)
        }

    q: dict[tuple[Any, Any], float] = {}
    alpha: dict[tuple[Any, Any], float] = {}
    for edge, prob in edge_probs.items():
        key = _canonical_edge_key(edge)
        q[key] = _prob_float(prob) if output == "float" else 0.0
        alpha[key] = _prob_first_order_coeff(prob)
    mean_actual = sum(q.values()) / len(q) if q else 0.0
    return _PreparedEdgeProbabilities(
        q=q,
        alpha=alpha,
        edge_failure_rate=None,
        mean_actual_edge_probability=mean_actual,
    )


def _tree_risk(
    tree: nx.Graph,
    root: Any,
    prepared: _PreparedEdgeProbabilities,
    *,
    total_weight: float,
    output: Literal["float", "poly"],
):
    if total_weight <= 0:
        raise ValueError("total graph weight is zero")
    if tree.number_of_nodes() == 0:
        return _zero(output)

    parent = {root: None}
    parent_edge: dict[Any, tuple[Any, Any]] = {}
    order = [root]
    stack = [root]
    while stack:
        node = stack.pop()
        for neighbor in tree.neighbors(node):
            if neighbor == parent.get(node):
                continue
            parent[neighbor] = node
            parent_edge[neighbor] = edge_key(node, neighbor)
            stack.append(neighbor)
            order.append(neighbor)

    if output == "float":
        live_prob = {root: 1.0}
        connected_load = float(tree.nodes[root].get("weight", 0.0))
        for node in order[1:]:
            edge = parent_edge[node]
            q_edge = _combined_failure_float(tree.edges[edge], prepared)
            live_prob[node] = live_prob[parent[node]] * (1.0 - q_edge)
            connected_load += float(tree.nodes[node].get("weight", 0.0)) * live_prob[node]
            connected_load += float(tree.edges[edge].get("edge_weight", 0.0)) * live_prob[node]
        return 1.0 - connected_load / total_weight

    path_a1 = {root: 0.0}
    path_a2 = {root: 0.0}
    c1 = 0.0
    c2 = 0.0
    for node in order[1:]:
        edge = parent_edge[node]
        alpha = _combined_failure_alpha(tree.edges[edge], prepared)
        prev = parent[node]
        path_a2[node] = path_a2[prev] + path_a1[prev] * alpha
        path_a1[node] = path_a1[prev] + alpha
        load = float(tree.nodes[node].get("weight", 0.0)) + float(tree.edges[edge].get("edge_weight", 0.0))
        c1 += load * path_a1[node]
        c2 -= load * path_a2[node]
    return PROBS.Poly([0.0, c1 / total_weight, c2 / total_weight])


def _nonbridge_section_edge_weight_risk(
    graph: nx.Graph,
    bridges: Iterable[tuple[Any, Any]],
    prepared: _PreparedEdgeProbabilities,
    *,
    total_weight: float,
    output: Literal["float", "poly"],
):
    if total_weight <= 0:
        raise ValueError("total graph weight is zero")
    bridge_set = {edge_key(*edge) for edge in bridges}
    c1 = 0.0
    for u, v, data in graph.edges(data=True):
        if edge_key(u, v) in bridge_set:
            continue
        edge_weight = float(data.get("edge_weight", 0.0))
        if edge_weight == 0.0:
            continue
        if output == "float":
            c1 += edge_weight * _combined_failure_float(data, prepared)
        else:
            c1 += edge_weight * _combined_failure_alpha(data, prepared)
    if output == "float":
        return c1 / total_weight
    return PROBS.Poly([0.0, c1 / total_weight])


def _chain_two_cut_risk(
    chains: Iterable[ChainSummary],
    graph: nx.Graph,
    prepared: _PreparedEdgeProbabilities,
    *,
    total_weight: float,
    output: Literal["float", "poly"],
):
    if total_weight <= 0:
        raise ValueError("total graph weight is zero")
    coefficient = 0.0
    for chain in chains:
        coefficient += _chain_two_cut_damage_sum(chain, graph, prepared, output=output)
    if output == "float":
        return coefficient / total_weight
    return PROBS.Poly([0.0, 0.0, coefficient / total_weight])


def _chain_two_cut_damage_sum(
    chain: ChainSummary,
    graph: nx.Graph,
    prepared: _PreparedEdgeProbabilities,
    *,
    output: Literal["float", "poly"],
) -> float:
    if len(chain.edges) < 2:
        return 0.0

    probs = [
        _combined_failure_float(graph.edges[edge], prepared)
        if output == "float"
        else _combined_failure_alpha(graph.edges[edge], prepared)
        for edge in chain.edges
    ]
    prefix_load = _chain_prefix_loads(chain, graph)
    left_prob_sum = 0.0
    left_prob_prefix_sum = 0.0
    damage = 0.0
    for idx, prob in enumerate(probs):
        damage += prob * (prefix_load[idx] * left_prob_sum - left_prob_prefix_sum)
        left_prob_sum += prob
        left_prob_prefix_sum += prob * prefix_load[idx]
    return damage


def _chain_prefix_loads(chain: ChainSummary, graph: nx.Graph) -> list[float]:
    # prefix[t] is load between the chain start and node position t:
    # nodes 1..t plus edge loads 1..t-1. Two failed edges i<j trap
    # prefix[j] - prefix[i].
    prefix = [0.0]
    acc = 0.0
    for pos in range(1, len(chain.nodes)):
        node = chain.nodes[pos]
        acc += float(graph.nodes[node].get("weight", 0.0))
        if pos - 1 > 0:
            acc += float(graph.edges[chain.edges[pos - 1]].get("edge_weight", 0.0))
        prefix.append(acc)
    return prefix


def _combined_failure_float(edge_data: dict, prepared: _PreparedEdgeProbabilities) -> float:
    probs = [_edge_q(edge, prepared) for edge in edge_data.get("original_edges", ())]
    if not probs:
        probs = [_edge_q(edge_data.get("original_edge"), prepared)] if edge_data.get("original_edge") is not None else []
    if not probs:
        return 0.0
    live = 1.0
    for prob in probs:
        live *= 1.0 - prob
    return 1.0 - live


def _combined_failure_alpha(edge_data: dict, prepared: _PreparedEdgeProbabilities) -> float:
    edges = edge_data.get("original_edges", ())
    if not edges and edge_data.get("original_edge") is not None:
        edges = (edge_data["original_edge"],)
    return sum(_edge_alpha(edge, prepared) for edge in edges)


def _edge_q(edge: Any, prepared: _PreparedEdgeProbabilities) -> float:
    return prepared.q.get(_canonical_edge_key(edge), 0.0)


def _edge_alpha(edge: Any, prepared: _PreparedEdgeProbabilities) -> float:
    return prepared.alpha.get(_canonical_edge_key(edge), 0.0)


def _canonical_edge_key(edge: Any) -> tuple[Any, Any]:
    if edge is None:
        return (None, None)
    if len(edge) >= 2:
        return edge_key(edge[0], edge[1])
    raise ValueError(f"edge key must contain at least two endpoints, got {edge!r}")


def _prob_float(prob: Any) -> float:
    if isinstance(prob, PROBS.Float):
        return float(prob)
    if isinstance(prob, PROBS.Poly):
        if prob.degree == 0:
            return float(prob[0])
        raise ValueError("cannot evaluate a non-constant Poly probability without an explicit p value")
    if hasattr(prob, "prob"):
        return float(prob.prob)
    return float(prob)


def _prob_first_order_coeff(prob: Any) -> float:
    if isinstance(prob, PROBS.Poly):
        return float(prob[1])
    if isinstance(prob, PROBS.Float):
        return float(prob)
    if hasattr(prob, "prob"):
        return float(prob.prob)
    return float(prob)


def _zero(output: Literal["float", "poly"]):
    return 0.0 if output == "float" else PROBS.Poly([0.0])


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


def _contract_three_edge_components(
    graph: nx.Graph,
    sources: Iterable[Any],
    *,
    chains: Iterable[ChainSummary] | None = None,
    component_method: ThreeEdgeComponentMethod = "networkx",
) -> nx.Graph:
    if graph.number_of_nodes() == 0:
        return nx.Graph()

    source_set = set(sources)
    components = _three_edge_components_from_chains(
        graph,
        chains=chains,
        component_method=component_method,
    )
    node_to_component: dict[Any, Any] = {}
    component_nodes: list[tuple[Any, set[Any]]] = []
    for idx, comp in enumerate(components):
        if len(comp) < 2:
            continue
        overlap = comp & set(node_to_component)
        if overlap:
            raise ValueError(
                "3-edge component detection returned overlapping components; "
                f"examples: {sorted(overlap, key=repr)[:5]}"
            )
        macro_node = _component_node_name("C3", idx, comp, source_set)
        component_nodes.append((macro_node, comp))
        for node in comp:
            node_to_component[node] = macro_node

    macro = nx.Graph()
    for macro_node, comp in component_nodes:
        original_nodes = _collect_original_nodes(graph, comp)
        internal_edge_weight = sum(
            float(data.get("edge_weight", 0.0))
            for u, v, data in graph.edges(comp, data=True)
            if u in comp and v in comp
        )
        macro.add_node(
            macro_node,
            weight=sum(float(graph.nodes[node].get("weight", 0.0)) for node in comp) + internal_edge_weight,
            original_nodes=frozenset(original_nodes),
            members=frozenset(comp),
            contains_source=bool(comp & source_set),
        )
    for node, data in graph.nodes(data=True):
        if node in node_to_component:
            continue
        node_data = data.copy()
        node_data.setdefault("original_nodes", frozenset({node}))
        node_data.setdefault("members", frozenset({node}))
        node_data["contains_source"] = node in source_set
        macro.add_node(node, **node_data)
        node_to_component[node] = node

    for u, v, data in graph.edges(data=True):
        cu = node_to_component[u]
        cv = node_to_component[v]
        if cu == cv:
            continue
        if macro.has_edge(cu, cv):
            _merge_macro_edge_data(macro.edges[cu, cv], data)
            continue
        edge_data = _initial_macro_edge_data(data)
        macro.add_edge(cu, cv, **edge_data)
    return macro


def _initial_macro_edge_data(data: dict) -> dict:
    edge_data = data.copy()
    edge_data.setdefault("edge_weight", 0.0)
    edge_data.setdefault("length", data.get("length_m", 1.0))
    edge_data.setdefault("parallel_macro_edge_count", 1)
    edge_data.setdefault("parallel_macro_edges", (tuple(edge_data.get("original_edges", ())),))
    return edge_data


def _merge_macro_edge_data(existing: dict, new_data: dict) -> None:
    existing["parallel_macro_edge_count"] = int(existing.get("parallel_macro_edge_count", 1)) + 1
    existing["length"] = float(existing.get("length", 0.0)) + float(new_data.get("length", new_data.get("length_m", 1.0)))
    existing["edge_weight"] = float(existing.get("edge_weight", 0.0)) + float(new_data.get("edge_weight", 0.0))
    existing["original_nodes"] = frozenset(
        set(existing.get("original_nodes", set())) | set(new_data.get("original_nodes", set()))
    )
    existing["original_edges"] = frozenset(
        set(existing.get("original_edges", set())) | set(new_data.get("original_edges", set()))
    )
    existing["parallel_macro_edges"] = tuple(existing.get("parallel_macro_edges", ())) + (
        tuple(new_data.get("original_edges", ())),
    )


def _three_edge_components_from_chains(
    graph: nx.Graph,
    *,
    chains: Iterable[ChainSummary] | None,
    component_method: ThreeEdgeComponentMethod,
) -> list[set[Any]]:
    if component_method not in ("networkx", "projection"):
        raise ValueError("component_method must be 'networkx' or 'projection'")

    analysis_graph = _analysis_graph_from_chains(graph, chains)
    if component_method == "networkx":
        return _raw_three_edge_components_networkx(analysis_graph)
    return _raw_three_edge_components_projection(analysis_graph)


def _analysis_graph_from_chains(
    graph: nx.Graph,
    chains: Iterable[ChainSummary] | None,
) -> nx.Graph:
    analysis = nx.Graph()
    analysis.add_nodes_from(graph.nodes)

    if chains is None:
        chain_iter = (
            ChainSummary(
                endpoints=(u, v),
                nodes=(u, v),
                edges=(edge_key(u, v),),
                internal_nodes=(),
                length=float(data.get("length", data.get("length_m", 1.0))),
                weight=0.0,
                edge_weight=float(data.get("edge_weight", 0.0)),
            )
            for u, v, data in graph.edges(data=True)
        )
    else:
        chain_iter = iter(chains)

    for chain_id, chain in enumerate(chain_iter):
        u, v = chain.endpoints
        analysis.add_node(u)
        analysis.add_node(v)
        if u == v:
            continue
        aux = ("reduced_edge", chain_id)
        analysis.add_edge(u, aux)
        analysis.add_edge(aux, v)
    return analysis


def _raw_three_edge_components_networkx(analysis_graph: nx.Graph) -> list[set[Any]]:
    raw_components: list[set[Any]] = []
    for component in nx.k_edge_components(analysis_graph, k=3):
        skeleton_nodes = {node for node in component if not isinstance(node, tuple)}
        if len(skeleton_nodes) >= 2:
            raw_components.append(set(skeleton_nodes))
    return raw_components


def _raw_three_edge_components_projection(analysis_graph: nx.Graph) -> list[set[Any]]:
    cut_detector = _load_cascading_projection_cut_detector()
    raw_components: list[set[Any]] = []
    for nodes in nx.connected_components(analysis_graph):
        subgraph = analysis_graph.subgraph(nodes).copy()
        skeleton_nodes = {node for node in subgraph.nodes if not _is_reduced_edge_node(node)}
        if len(skeleton_nodes) < 2:
            continue
        raw_components.extend(
            _raw_three_edge_components_from_two_cuts(
                subgraph,
                cut_detector(subgraph),
                min_skeleton_nodes=2,
            )
        )
    return raw_components


def _raw_three_edge_components_from_two_cuts(
    analysis_graph: nx.Graph,
    two_edge_cuts: Iterable[tuple[tuple[Any, Any], tuple[Any, Any]]],
    *,
    min_skeleton_nodes: int,
) -> list[set[Any]]:
    """Return skeleton 3-edge components from algebraic two-edge-cut classes.

    In a bridgeless component, two-edge cuts form edge equivalence classes.
    Removing all nontrivial cut-class edges leaves the 3-edge-connected
    skeleton blocks. This avoids copying the full graph once per detected cut.
    """
    edge_parent = {edge_key(u, v): edge_key(u, v) for u, v in analysis_graph.edges}

    def find(edge: tuple[Any, Any]) -> tuple[Any, Any]:
        parent = edge_parent[edge]
        if parent != edge:
            edge_parent[edge] = find(parent)
        return edge_parent[edge]

    def union(edge_a: tuple[Any, Any], edge_b: tuple[Any, Any]) -> None:
        root_a = find(edge_a)
        root_b = find(edge_b)
        if root_a != root_b:
            edge_parent[root_b] = root_a

    for raw_a, raw_b in two_edge_cuts:
        if _share_reduced_edge_node(raw_a, raw_b):
            continue
        edge_a = edge_key(*raw_a)
        edge_b = edge_key(*raw_b)
        if edge_a not in edge_parent or edge_b not in edge_parent:
            continue
        union(edge_a, edge_b)

    classes: dict[tuple[Any, Any], list[tuple[Any, Any]]] = {}
    for edge in edge_parent:
        classes.setdefault(find(edge), []).append(edge)
    separator_edges = [
        edge
        for edges in classes.values()
        if len(edges) >= 2
        for edge in edges
    ]

    reduced = analysis_graph.copy()
    reduced.remove_edges_from(separator_edges)
    components: list[set[Any]] = []
    for nodes in nx.connected_components(reduced):
        skeleton = {node for node in nodes if not _is_reduced_edge_node(node)}
        if len(skeleton) >= min_skeleton_nodes:
            components.append(skeleton)
    return components


def _is_reduced_edge_node(node: Any) -> bool:
    return isinstance(node, tuple) and len(node) >= 1 and node[0] == "reduced_edge"


def _share_reduced_edge_node(edge_a: tuple[Any, Any], edge_b: tuple[Any, Any]) -> bool:
    return bool(
        {_node for _node in edge_a if _is_reduced_edge_node(_node)}
        & {_node for _node in edge_b if _is_reduced_edge_node(_node)}
    )


def _load_cascading_projection_cut_detector():
    try:
        from dc_graph.structure import _near_zero_two_edge_cuts_from_projection

        return _near_zero_two_edge_cuts_from_projection
    except ModuleNotFoundError:
        sibling = Path(__file__).resolve().parents[2] / "cascading"
        if sibling.exists():
            sys.path.insert(0, str(sibling))
            from dc_graph.structure import _near_zero_two_edge_cuts_from_projection

            return _near_zero_two_edge_cuts_from_projection
        raise ModuleNotFoundError(
            "Projection-based 3-edge decomposition requires the sibling "
            f"cascading project at {sibling}"
        )


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
