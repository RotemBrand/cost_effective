import networkx as nx

from indexes.switch_placement import (
    add_synthetic_switches,
    count_switch_edges,
    count_tie_edges,
)


def _square_with_diagonal():
    graph = nx.Graph()
    graph.add_node(0, pos=(0.0, 0.0))
    graph.add_node(1, pos=(1.0, 0.0))
    graph.add_node(2, pos=(1.0, 1.0))
    graph.add_node(3, pos=(0.0, 1.0))
    graph.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])
    return graph


def test_add_synthetic_switches_keeps_topology_and_marks_non_tree_ties():
    graph = _square_with_diagonal()

    switched = add_synthetic_switches(graph, n_switches=3, seed=1)

    assert set(switched.edges) == set(graph.edges)
    assert switched.number_of_edges() - switched.number_of_nodes() + 1 == 2
    assert count_tie_edges(switched) == 2
    assert count_switch_edges(switched) == 3


def test_add_synthetic_switches_keeps_all_ties_even_if_switch_budget_is_smaller():
    graph = _square_with_diagonal()

    switched = add_synthetic_switches(graph, n_switches=1, seed=1)

    assert count_tie_edges(switched) == 2
    assert count_switch_edges(switched) == 2


def test_add_synthetic_switches_is_seed_reproducible():
    graph = nx.path_graph(8)
    for node in graph:
        graph.nodes[node]["pos"] = (float(node), 0.0)
    graph.add_edge(0, 7)

    switched_a = add_synthetic_switches(graph, n_switches=4, seed=10)
    switched_b = add_synthetic_switches(graph, n_switches=4, seed=10)

    attrs_a = nx.get_edge_attributes(switched_a, "is_switch")
    attrs_b = nx.get_edge_attributes(switched_b, "is_switch")
    assert attrs_a == attrs_b
