import networkx as nx
import pytest

from indexes import Float, GraphRel, contract_switch_sections
from indexes.simulation import ConnectedComponents


def _set_node_weights(graph, weights):
    nx.set_node_attributes(graph, weights, "weight")
    return graph


def test_zero_edge_weight_matches_omitted_edge_weight_mcmc():
    graph = nx.path_graph(4)
    _set_node_weights(graph, {0: 0.0, 1: 1.0, 2: 2.0, 3: 3.0})
    edge_probs = {edge: Float(0.05) for edge in graph.edges}

    without_edge_weight = GraphRel(graph, edges_prob=edge_probs, sources=[0])
    with_zero_edge_weight = GraphRel(graph, edges_prob=edge_probs, edges_weight=0.0, sources=[0])

    res_without = without_edge_weight.calc_rel_simulation(
        T_days=30,
        mean_cycle_days=0.5,
        seed=4,
    )
    res_with_zero = with_zero_edge_weight.calc_rel_simulation(
        T_days=30,
        mean_cycle_days=0.5,
        seed=4,
    )

    assert res_with_zero.rel_result == pytest.approx(res_without.rel_result)


def test_failed_edge_section_weight_is_disconnected():
    graph = nx.Graph()
    graph.add_edge(0, 1, edge_weight=5.0)
    _set_node_weights(graph, {0: 0.0, 1: 0.0})

    components = ConnectedComponents(
        graph,
        source=0,
        weight_attr="weight",
        rel_type="saidi",
        edge_weight_attr="edge_weight",
    )

    assert components.reliability() == pytest.approx(0.0)
    components.remove_edge((0, 1))
    assert components.reliability() == pytest.approx(1.0)


def test_downstream_live_edge_weight_is_disconnected_after_upstream_failure():
    graph = nx.Graph()
    graph.add_edge(0, 1, edge_weight=2.0)
    graph.add_edge(1, 2, edge_weight=3.0)
    _set_node_weights(graph, {0: 0.0, 1: 0.0, 2: 0.0})

    components = ConnectedComponents(
        graph,
        source=0,
        weight_attr="weight",
        rel_type="saidi",
        edge_weight_attr="edge_weight",
    )

    components.remove_edge((0, 1))
    assert components.reliability() == pytest.approx(1.0)

    components = ConnectedComponents(
        graph,
        source=0,
        weight_attr="weight",
        rel_type="saidi",
        edge_weight_attr="edge_weight",
    )
    components.remove_edge((1, 2))
    assert components.reliability() == pytest.approx(3.0 / 5.0)


def test_parallel_tie_preserves_connectivity_but_not_failed_section_load():
    graph = nx.MultiGraph()
    graph.add_edge(0, 1, key="section", edge_weight=5.0, is_tie=False)
    graph.add_edge(0, 1, key="tie", edge_weight=0.0, is_tie=True)
    _set_node_weights(graph, {0: 0.0, 1: 2.0})

    components = ConnectedComponents(
        graph,
        source=0,
        weight_attr="weight",
        rel_type="saidi",
        edge_weight_attr="edge_weight",
    )

    components.remove_edge((0, 1, "section"))
    assert components.reliability() == pytest.approx(5.0 / 7.0)

    components = ConnectedComponents(
        graph,
        source=0,
        weight_attr="weight",
        rel_type="saidi",
        edge_weight_attr="edge_weight",
    )
    components.remove_edge((0, 1, "tie"))
    assert components.reliability() == pytest.approx(0.0)


def test_contract_switch_sections_keeps_tie_parallel_to_contracted_section():
    graph = nx.Graph()
    graph.add_edge(0, 1, length=1.0, prob=0.1)
    graph.add_edge(1, 2, length=2.0, prob=0.2)
    graph.add_edge(2, 3, length=3.0, prob=0.3)
    graph.add_edge(1, 3, is_tie=True, length=4.0, prob=0.4)
    _set_node_weights(graph, {0: 0.0, 1: 0.0, 2: 7.0, 3: 0.0})

    section_graph = contract_switch_sections(graph, sources=[0])

    assert isinstance(section_graph, nx.MultiGraph)
    assert 2 not in section_graph
    assert section_graph.number_of_edges(1, 3) == 2

    edges_1_3 = list(section_graph.get_edge_data(1, 3).values())
    section_edges = [data for data in edges_1_3 if not data["is_tie"]]
    tie_edges = [data for data in edges_1_3 if data["is_tie"]]

    assert len(section_edges) == 1
    assert len(tie_edges) == 1
    assert section_edges[0]["edge_weight"] == pytest.approx(7.0)
    assert section_edges[0]["length"] == pytest.approx(5.0)
    assert section_edges[0]["prob"] == pytest.approx(1 - (1 - 0.2) * (1 - 0.3))


def test_contract_switch_sections_uses_closed_switches_as_boundaries():
    graph = nx.path_graph(5)
    nx.set_node_attributes(graph, 1.0, "weight")
    nx.set_edge_attributes(graph, 0.1, "prob")
    nx.set_edge_attributes(graph, 1.0, "length")
    nx.set_edge_attributes(graph, False, "is_tie")
    nx.set_edge_attributes(graph, False, "is_switch")
    graph.edges[2, 3]["is_switch"] = True

    section_graph = contract_switch_sections(graph, sources=[0])

    assert set(section_graph.nodes) == {0, 2, 3, 4}
    original_edges = [
        data["original_edges"]
        for _, _, data in section_graph.edges(data=True)
        if not data["is_tie"]
    ]
    assert [(0, 1), (1, 2)] in original_edges
    assert [(2, 3)] in original_edges
    assert [(3, 4)] in original_edges


def test_exact_reliability_rejects_nonzero_edge_weights():
    graph = nx.path_graph(3)
    _set_node_weights(graph, {0: 0.0, 1: 1.0, 2: 1.0})
    graph.edges[0, 1]["edge_weight"] = 1.0

    graph_rel = GraphRel(graph, sources=[0])

    with pytest.raises(NotImplementedError):
        graph_rel.calc_rel(["saidi"])
