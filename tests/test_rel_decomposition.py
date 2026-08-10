import networkx as nx
import numpy as np
import pytest

from indexes import GraphRel
from indexes.rel_decomposition import extract_chains
from indexes.probs import Float, Poly


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


def _two_blobs_with_two_macro_links() -> nx.Graph:
    graph = nx.Graph()
    graph.add_edges_from(nx.complete_graph([0, 1, 2, 3]).edges())
    graph.add_edges_from(nx.complete_graph([4, 5, 6, 7]).edges())
    graph.add_edges_from([(2, 4), (3, 5)])
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


def test_generalized_macro_aggregates_parallel_links_in_simple_graph():
    decomposition = _graph_rel(_two_blobs_with_two_macro_links(), sources=[0]).decompose()
    parallel_counts = [
        int(data.get("parallel_macro_edge_count", 1))
        for _, _, data in decomposition.three_edge_macro_graph.edges(data=True)
    ]

    assert decomposition.three_edge_macro_graph.number_of_edges() == 1
    assert parallel_counts == [2]


def test_projection_generalized_components_match_networkx():
    graph = _three_blob_cycle_graph()
    networkx_decomp = _graph_rel(graph, sources=[0]).decompose(generalized_component_method="networkx")
    projection_decomp = _graph_rel(graph, sources=[0]).decompose(generalized_component_method="projection")

    networkx_members = {
        frozenset(data["members"])
        for _, data in networkx_decomp.three_edge_macro_graph.nodes(data=True)
    }
    projection_members = {
        frozenset(data["members"])
        for _, data in projection_decomp.three_edge_macro_graph.nodes(data=True)
    }

    assert projection_members == networkx_members
    assert len(projection_decomp.generalized_chains) == len(networkx_decomp.generalized_chains)


def test_generalized_chain_structural_risk_uses_three_edge_macro_graph():
    graph = _three_blob_cycle_graph()
    probs = {edge: Poly([0, 1]) for edge in graph.edges}

    decomp = _graph_rel(graph, sources=[0], edges_prob=probs).decompose(generalized_component_method="projection")
    terms = decomp.switch_risk_terms(output="poly")

    assert decomp.three_edge_macro_graph.number_of_nodes() == 3
    assert decomp.source_graph.number_of_nodes() == 12
    assert terms.structural[2] == pytest.approx(16 / 11)


def test_tree_td_matches_full_td_on_pure_tree():
    graph = nx.path_graph(4)
    _set_weights(graph, {0: 0, 1: 2, 2: 3, 3: 4})
    decomposition = _graph_rel(graph, sources=[0], max_fail=2).decompose()

    full = decomposition.td_saidi(max_fail=2).prob
    tree = decomposition.tree_td_saidi(max_fail=2).prob

    assert np.allclose(full, tree)


def test_switch_tree_risk_includes_bridge_edge_weights_float_and_poly():
    graph = nx.path_graph(3)
    _set_weights(graph, {0: 0, 1: 0, 2: 10})
    graph.edges[0, 1]["edge_weight"] = 3.0
    graph.edges[1, 2]["edge_weight"] = 4.0
    float_probs = {(0, 1): Float(0.1), (1, 2): Float(0.2)}
    poly_probs = {(0, 1): Poly([0, 1]), (1, 2): Poly([0, 2])}

    decomp = _graph_rel(graph, sources=[0], edges_prob=float_probs).decompose(include_generalized_chains=False)
    float_terms = decomp.switch_risk_terms(output="float")
    poly_terms = decomp.switch_risk_terms(edge_probs=poly_probs, output="poly")

    assert decomp.total_weight == pytest.approx(17.0)
    assert float_terms.tree == pytest.approx((3 * 0.1 + 14 * (1 - 0.9 * 0.8)) / 17)
    assert poly_terms.tree[1] == pytest.approx(45 / 17)
    assert poly_terms.tree[2] == pytest.approx(-28 / 17)


def test_regular_chain_two_cut_risk_is_exact_p2_coefficient():
    graph = nx.cycle_graph(4)
    _set_weights(graph, {0: 0, 1: 2, 2: 3, 3: 5})
    probs = {edge: Poly([0, 1]) for edge in graph.edges}

    decomp = _graph_rel(graph, sources=[0], edges_prob=probs).decompose(include_generalized_chains=False)
    terms = decomp.switch_risk_terms(output="poly")

    assert len(decomp.regular_chains) == 1
    assert terms.internal[1] == pytest.approx(0.0)
    assert terms.internal[2] == pytest.approx(33 / 10)


def test_length_mean_failure_probability_is_calibrated():
    graph = nx.cycle_graph(4)
    _set_weights(graph, {0: 0})
    for idx, edge in enumerate(graph.edges, start=1):
        graph.edges[edge]["length"] = float(idx)

    decomp = _graph_rel(graph, sources=[0]).decompose(include_generalized_chains=False)
    terms = decomp.switch_risk_terms(mean_edge_failure_prob=5e-4, length_attr="length", output="float")

    assert terms.mean_actual_edge_probability == pytest.approx(5e-4, abs=1e-5)
    assert terms.edge_failure_rate is not None
