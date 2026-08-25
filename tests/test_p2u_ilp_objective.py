import networkx as nx
import pytest

gp = pytest.importorskip("gurobipy")

from figures.optimal_sfo.run_p2u_ilp_2edge import solve_min_2edge


def _complete_four_with_one_loaded_section() -> nx.Graph:
    graph = nx.Graph()
    for node in range(4):
        graph.add_node(str(node), size_kva=1.0)

    edge_data = {
        ("0", "1"): (1.0, 0.0),
        ("1", "2"): (1.0, 0.0),
        ("2", "3"): (1.0, 0.0),
        ("0", "3"): (1.0, 0.0),
        ("0", "2"): (2.0, 0.0),
        ("1", "3"): (100.0, 10.0),
    }
    for idx, ((u, v), (length, demand)) in enumerate(edge_data.items()):
        graph.add_edge(
            u,
            v,
            edge_id=idx,
            length_m=length,
            edge_size_kva=demand,
            edge_transformer_count=int(demand > 0),
        )
    return graph


def test_chain_demand_objective_can_choose_loaded_contracted_edge():
    graph = _complete_four_with_one_loaded_section()

    min_length, min_summary = solve_min_2edge(
        graph,
        original_source_incident=None,
        redundancy=2,
        max_redundancy=None,
        time_limit=10.0,
        mip_gap=0.0,
        threads=1,
        max_cut_rounds=20,
        cut_mode="iterative",
        objective_mode="min_length",
    )
    demand_first, demand_summary = solve_min_2edge(
        graph,
        original_source_incident=None,
        redundancy=2,
        max_redundancy=None,
        time_limit=10.0,
        mip_gap=0.0,
        threads=1,
        max_cut_rounds=20,
        cut_mode="iterative",
        objective_mode="max_chain_demand_then_min_length",
    )

    assert min_length is not None
    assert demand_first is not None
    assert min_summary["selected_edge_coverage_weight"] == pytest.approx(0.0)
    assert demand_summary["selected_edge_coverage_weight"] == pytest.approx(10.0)
    assert demand_first.has_edge("1", "3")
    assert demand_summary["selected_edge_coverage_fraction"] == pytest.approx(1.0)
    assert demand_summary["objective_length_m"] > min_summary["objective_length_m"]
