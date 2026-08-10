import networkx as nx
import numpy as np
import pytest

from indexes import GraphRel
from indexes.rel_decomposition import extract_chains


def _set_weights(graph: nx.Graph, weights: dict | None = None) -> nx.Graph:
    weights = weights or {}
    for node in graph.nodes:
        graph.nodes[node]["weight"] = float(weights.get(node, 1.0))
    for edge in graph.edges:
        graph.edges[edge]["length"] = 1.0
        graph.edges[edge]["edge_weight"] = 0.0
    return graph


def _graph_rel(graph: nx.Graph, **kwargs) -> GraphRel:
    return GraphRel(graph, nodes_weight=nx.get_node_attributes(graph, "weight"), **kwargs)


def test_decompose_rejects_multigraph_input():
    graph = nx.MultiGraph()
    graph.add_edge(0, 1)
    graph.add_edge(0, 1)
    _set_weights(graph)

    graph_rel = _graph_rel(graph, sources=[0])

    with pytest.raises(TypeError, match="nx.Graph input only"):
        graph_rel.decompose()


def test_decompose_rejects_parallel_edges_created_by_source_contraction():
    graph = nx.Graph()
    graph.add_edges_from([(0, 2), (1, 2), (2, 3)])
    _set_weights(graph)

    graph_rel = _graph_rel(graph, sources=[0, 1])

    with pytest.raises(TypeError, match="parallel edges after source contraction"):
        graph_rel.decompose()


def test_bridge_tree_aggregates_two_edge_block_weights():
    graph = nx.Graph()
    graph.add_edges_from([(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)])
    _set_weights(graph, {0: 0, 1: 2, 2: 3, 3: 5, 4: 7})

    decomposition = _graph_rel(graph, sources=[0]).decompose()

    bridge_weights = sorted(data["weight"] for _, data in decomposition.bridge_tree.nodes(data=True))
    assert bridge_weights == [5.0, 5.0, 7.0]
    assert len(decomposition.bridges) == 2


def test_extract_chains_preserves_cycle_internal_weight():
    graph = nx.cycle_graph(4)
    _set_weights(graph, {0: 0, 1: 2, 2: 3, 3: 4})

    chains = extract_chains(graph, sources=[0])

    assert len(chains) == 1
    chain = chains[0]
    assert chain.endpoints == (0, 0)
    assert set(chain.internal_nodes) == {1, 2, 3}
    assert chain.weight == pytest.approx(9.0)
    assert chain.length == pytest.approx(4.0)


def _three_blob_cycle_graph() -> nx.Graph:
    graph = nx.Graph()
    graph.add_edges_from(nx.complete_graph([0, 1, 2, 3]).edges())
    graph.add_edges_from(nx.complete_graph([4, 5, 6, 7]).edges())
    graph.add_edges_from(nx.complete_graph([8, 9, 10, 11]).edges())
    graph.add_edges_from([(2, 4), (6, 8), (10, 3)])
    _set_weights(graph, {0: 0})
    return graph


def test_generalized_chains_contract_three_edge_components():
    decomposition = _graph_rel(_three_blob_cycle_graph(), sources=[0]).decompose()

    assert decomposition.three_edge_macro_graph.number_of_nodes() == 3
    assert decomposition.three_edge_macro_graph.number_of_edges() == 3
    assert len(decomposition.generalized_chains) == 1
    generalized = decomposition.generalized_chains[0]
    assert generalized.endpoints == (0, 0)
    assert len(generalized.internal_nodes) == 2
    assert generalized.weight == pytest.approx(8.0)


def test_tree_td_matches_full_td_on_pure_tree():
    graph = nx.path_graph(4)
    _set_weights(graph, {0: 0, 1: 2, 2: 3, 3: 4})
    decomposition = _graph_rel(graph, sources=[0], max_fail=2).decompose()

    full = decomposition.td_saidi(max_fail=2).prob
    tree = decomposition.tree_td_saidi(max_fail=2).prob

    assert np.allclose(full, tree)
