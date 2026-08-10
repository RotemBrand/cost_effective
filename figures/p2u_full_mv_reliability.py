from __future__ import annotations

import argparse
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "power" / "better_grids" / "SFO" / "P2U"
OUT_DIR = ROOT / "outputs" / "p2u_full_mv_reliability"
SOURCE_CRS = "EPSG:32610"


def _clean_node(node) -> str:
    return str(node).strip()


def _edge_key(u, v) -> tuple[str, str]:
    return tuple(sorted((_clean_node(u), _clean_node(v))))


def _read_layer(name: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(DATA_DIR / f"{name}.shp")
    if gdf.crs is None:
        gdf = gdf.set_crs(SOURCE_CRS)
    return gdf


def _substation_node_from_label(label: str) -> str:
    return _clean_node(label)


def build_full_mv_graph() -> tuple[nx.Graph, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Build the full P2U physical MV graph from Line_N.

    The graph contains normally closed lines and normally open ties. Edge
    `is_tie=True` means `Status == 0`.
    """
    lines = _read_layer("Line_N")
    mv = lines[lines["NomV"].astype(float).between(1, 40)].copy()
    switches = _read_layer("SwitchingDevices_N")
    switches = switches[switches["NomV_kV"].astype(float).between(1, 40)].copy()
    switch_edges = {_edge_key(row.NodeA, row.NodeB) for _, row in switches.iterrows()}

    graph = nx.Graph()
    for _, row in mv.iterrows():
        u = _clean_node(row.NodeA)
        v = _clean_node(row.NodeB)
        status = str(row.Status).strip()
        is_tie = status == "0"
        substation = _clean_node(row.Subest)
        length_m = float(row.geometry.length)
        if graph.has_edge(u, v):
            # Rare duplicate line features exist. Keep a single reliability edge
            # with summed length and conservative metadata.
            data = graph.edges[u, v]
            data["length_m"] += length_m
            data["raw_line_count"] = int(data.get("raw_line_count", 1)) + 1
            data["is_tie"] = data["is_tie"] or is_tie
            data["is_switch"] = data["is_switch"] or is_tie or (_edge_key(u, v) in switch_edges)
            data["statuses"] = sorted(set(data["statuses"]) | {status})
            continue

        graph.add_edge(
            u,
            v,
            length_m=length_m,
            length=length_m,
            raw_len_km=float(row.Len),
            is_tie=is_tie,
            is_switch=is_tie or (_edge_key(u, v) in switch_edges),
            normally_closed=not is_tie,
            status=int(status) if status.isdigit() else status,
            statuses=[status],
            substation=None if substation in {"", "True", "nan", "None"} else substation,
            feeder=None if _clean_node(row.Feeder) in {"", "True", "nan", "None"} else _clean_node(row.Feeder),
            consumer_count=0.0,
            demand_kw=0.0,
            yearly_kwh=0.0,
            is_source=False,
            raw_line_count=1,
        )
        for node in (u, v):
            graph.nodes[node].setdefault("weight", 0.0)
            graph.nodes[node].setdefault("consumer_count", 0.0)
            graph.nodes[node].setdefault("demand_kw", 0.0)
            graph.nodes[node].setdefault("yearly_kwh", 0.0)
            graph.nodes[node].setdefault("is_source", False)

    hvmv = _read_layer("HVMVSubstation_N")
    source_nodes: list[str] = []
    for _, row in hvmv.iterrows():
        source = f"{_clean_node(row.Node)}_1247"
        if source in graph:
            graph.nodes[source]["is_source"] = True
            source_nodes.append(source)
    graph.graph["sources"] = source_nodes
    return graph, mv, switches


def aggregate_lv_consumers_to_transformers(graph: nx.Graph) -> dict:
    """Aggregate LV consumer demand to MV transformer nodes using LV components."""
    lv_lines = _read_layer("Line_N")
    lv_lines = lv_lines[lv_lines["NomV"].astype(float).between(0.0, 1.0)].copy()
    consumers = _read_layer("NewConsumerGreenfield_N")
    transformers = _read_layer("DistribTransf_N")

    lv_graph = nx.Graph()
    for _, row in lv_lines.iterrows():
        lv_graph.add_edge(_clean_node(row.NodeA), _clean_node(row.NodeB))

    transformer_nodes = {_clean_node(node) for node in transformers["Node"]}
    transformer_lv_to_mv = {f"{node}LV": node for node in transformer_nodes}
    transformers = transformers.reset_index(drop=True)
    transformer_coords = np.column_stack((transformers.geometry.x.values, transformers.geometry.y.values))
    transformer_tree = cKDTree(transformer_coords)
    transformer_node_by_idx = transformers["Node"].map(_clean_node).to_numpy()

    consumer_to_transformer: dict[str, str] = {}
    components_without_transformer = 0
    components_with_multiple_transformers = 0

    for comp in nx.connected_components(lv_graph):
        transformer_hits = sorted(
            transformer_lv_to_mv[node]
            for node in comp
            if node in transformer_lv_to_mv
        )
        if len(transformer_hits) == 0:
            components_without_transformer += 1
            continue
        if len(transformer_hits) > 1:
            components_with_multiple_transformers += 1
        transformer = transformer_hits[0]
        for node in comp:
            consumer_to_transformer[node] = transformer

    assigned_consumers = 0
    unassigned_consumers = 0
    fallback_nearest_transformer_consumers = 0
    fallback_distances_m: list[float] = []
    total_demand_kw = 0.0
    total_num_customers = 0.0
    total_yearly_kwh = 0.0

    transformer_stats = defaultdict(lambda: {"demand_kw": 0.0, "consumer_count": 0.0, "yearly_kwh": 0.0, "points": 0})
    for _, row in consumers.iterrows():
        consumer = _clean_node(row.Code)
        transformer = consumer_to_transformer.get(consumer)
        demand_kw = float(row.get("DemP_kW", 0.0) or 0.0)
        num_customers = float(row.get("NumCust", 0.0) or 0.0)
        yearly_kwh = float(row.get("Yearly_kWh", 0.0) or 0.0)
        total_demand_kw += demand_kw
        total_num_customers += num_customers
        total_yearly_kwh += yearly_kwh
        if transformer is None or transformer not in graph:
            distance, idx = transformer_tree.query([row.geometry.x, row.geometry.y], k=1)
            transformer = str(transformer_node_by_idx[int(idx)])
            fallback_nearest_transformer_consumers += 1
            fallback_distances_m.append(float(distance))
        if transformer not in graph:
            unassigned_consumers += 1
            continue
        assigned_consumers += 1
        transformer_stats[transformer]["demand_kw"] += demand_kw
        transformer_stats[transformer]["consumer_count"] += num_customers
        transformer_stats[transformer]["yearly_kwh"] += yearly_kwh
        transformer_stats[transformer]["points"] += 1

    for node, stats in transformer_stats.items():
        graph.nodes[node]["weight"] = stats["demand_kw"]
        graph.nodes[node]["demand_kw"] = stats["demand_kw"]
        graph.nodes[node]["consumer_count"] = stats["consumer_count"]
        graph.nodes[node]["yearly_kwh"] = stats["yearly_kwh"]
        graph.nodes[node]["consumer_points"] = stats["points"]

    return {
        "lv_components": nx.number_connected_components(lv_graph),
        "lv_nodes": lv_graph.number_of_nodes(),
        "lv_edges": lv_graph.number_of_edges(),
        "lv_components_without_transformer": components_without_transformer,
        "lv_components_with_multiple_transformers": components_with_multiple_transformers,
        "consumer_points_total": int(len(consumers)),
        "consumer_points_assigned_to_mv": int(assigned_consumers),
        "consumer_points_unassigned": int(unassigned_consumers),
        "consumer_points_assigned_by_nearest_transformer_fallback": int(fallback_nearest_transformer_consumers),
        "nearest_transformer_fallback_mean_distance_m": float(np.mean(fallback_distances_m)) if fallback_distances_m else 0.0,
        "nearest_transformer_fallback_max_distance_m": float(np.max(fallback_distances_m)) if fallback_distances_m else 0.0,
        "demand_kw_total_raw": total_demand_kw,
        "demand_kw_assigned_to_mv": sum(stats["demand_kw"] for stats in transformer_stats.values()),
        "num_customers_total_raw": total_num_customers,
        "num_customers_assigned_to_mv": sum(stats["consumer_count"] for stats in transformer_stats.values()),
        "yearly_kwh_total_raw": total_yearly_kwh,
        "yearly_kwh_assigned_to_mv": sum(stats["yearly_kwh"] for stats in transformer_stats.values()),
        "transformer_nodes_with_assigned_consumers": len(transformer_stats),
    }


def _contract_sources(graph: nx.Graph, sources: Iterable[str], source_label: str = "__SOURCE__") -> nx.Graph:
    sources = set(sources)
    mapping = {source: source_label for source in sources}
    contracted = nx.relabel_nodes(graph, mapping, copy=True)
    if contracted.has_node(source_label):
        contracted.nodes[source_label]["weight"] = 0.0
    contracted.remove_edges_from(nx.selfloop_edges(contracted))
    return contracted


def two_edge_backbone_metrics(graph: nx.Graph, sources: list[str]) -> tuple[set[str], dict]:
    """Compute source-referenced and global 2-edge-connected backbone metrics."""
    contracted = _contract_sources(graph, sources)
    source_label = "__SOURCE__"
    bridges = set(nx.bridges(contracted))
    bridge_keys = {_edge_key(u, v) for u, v in bridges}

    without_bridges = contracted.copy()
    without_bridges.remove_edges_from(bridges)
    components = list(nx.connected_components(without_bridges))
    source_component = next((set(comp) for comp in components if source_label in comp), {source_label})
    source_backbone = set(source_component) - {source_label}

    nontrivial_components = [set(comp) for comp in components if len(comp - {source_label}) > 1]
    any_cycle_nodes = set().union(*(comp - {source_label} for comp in nontrivial_components)) if nontrivial_components else set()
    largest_component_size = max((len(comp - {source_label}) for comp in nontrivial_components), default=0)

    return source_backbone, {
        "bridge_edges_after_source_contraction": len(bridge_keys),
        "source_2edge_backbone_nodes": len(source_backbone),
        "nodes_in_any_2edge_component": len(any_cycle_nodes),
        "largest_2edge_component_nodes": largest_component_size,
        "nontrivial_2edge_components": len(nontrivial_components),
    }


def consumer_backbone_metrics(graph: nx.Graph, source_backbone: set[str]) -> dict:
    load_nodes = [node for node, data in graph.nodes(data=True) if float(data.get("demand_kw", 0.0)) > 0]
    load_nodes_backbone = [node for node in load_nodes if node in source_backbone]

    total_transformer_load = sum(float(graph.nodes[node].get("demand_kw", 0.0)) for node in load_nodes)
    backbone_transformer_load = sum(float(graph.nodes[node].get("demand_kw", 0.0)) for node in load_nodes_backbone)
    total_customers = sum(float(graph.nodes[node].get("consumer_count", 0.0)) for node in load_nodes)
    backbone_customers = sum(float(graph.nodes[node].get("consumer_count", 0.0)) for node in load_nodes_backbone)

    return {
        "transformer_load_nodes": len(load_nodes),
        "transformer_load_nodes_2edge_to_source": len(load_nodes_backbone),
        "percent_transformer_load_nodes_2edge_to_source": _percent(len(load_nodes_backbone), len(load_nodes)),
        "demand_kw_total": total_transformer_load,
        "demand_kw_2edge_to_source": backbone_transformer_load,
        "percent_demand_kw_2edge_to_source": _percent(backbone_transformer_load, total_transformer_load),
        "num_customers_total": total_customers,
        "num_customers_2edge_to_source": backbone_customers,
        "percent_customers_2edge_to_source": _percent(backbone_customers, total_customers),
    }


def _percent(part: float, whole: float) -> float:
    if whole == 0:
        return 0.0
    return 100.0 * float(part) / float(whole)


def topology_metrics(graph: nx.Graph, mv_rows: gpd.GeoDataFrame, switches: gpd.GeoDataFrame) -> dict:
    closed_edges = [(u, v) for u, v, data in graph.edges(data=True) if not data.get("is_tie", False)]
    open_ties = [(u, v) for u, v, data in graph.edges(data=True) if data.get("is_tie", False)]
    closed_graph = graph.edge_subgraph(closed_edges).copy()
    physical_cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + nx.number_connected_components(graph)
    closed_cycle_rank = closed_graph.number_of_edges() - closed_graph.number_of_nodes() + nx.number_connected_components(closed_graph)
    return {
        "raw_mv_line_features": int(len(mv_rows)),
        "nodes": graph.number_of_nodes(),
        "lines_unique_edges": graph.number_of_edges(),
        "normally_closed_edges": len(closed_edges),
        "open_tie_edges_R": len(open_ties),
        "raw_open_tie_features": int((mv_rows["Status"].astype(str).str.strip() == "0").sum()),
        "mv_switching_device_features": int(len(switches)),
        "total_ties_cycle_rank_physical_graph": int(physical_cycle_rank),
        "closed_graph_cycle_rank": int(closed_cycle_rank),
        "physical_components": nx.number_connected_components(graph),
        "closed_operating_components": nx.number_connected_components(closed_graph),
    }


def stationary_mc_saidi(
    graph: nx.Graph,
    sources: list[str],
    *,
    p: float,
    n_samples: int,
    seed: int,
    show_progress_every: int = 0,
) -> dict:
    """Estimate stationary SAIDI as expected disconnected demand fraction.

    Each real MV edge is independently failed with probability `p`.
    Source nodes are contracted before sampling.
    """
    if not 0 <= p < 1:
        raise ValueError("p must satisfy 0 <= p < 1")
    if n_samples <= 0:
        raise ValueError("n_samples must be positive")

    rng = np.random.default_rng(seed)
    working = _contract_sources(graph, sources)
    source = "__SOURCE__"
    edges = list(working.edges())
    edge_data = {(u, v): working.edges[u, v].copy() for u, v in edges}
    real_edges = [(u, v) for u, v in edges if source not in (u, v)]
    m = len(real_edges)
    weights = {node: float(data.get("weight", 0.0)) for node, data in working.nodes(data=True)}
    total_weight = sum(weights.values())
    if total_weight <= 0:
        raise ValueError("total graph weight is zero; cannot compute SAIDI")

    disconnected = np.empty(n_samples, dtype=float)
    failed_edge_counts = np.empty(n_samples, dtype=np.int32)
    t0 = time.perf_counter()

    for i in range(n_samples):
        k = int(rng.binomial(m, p))
        failed_edge_counts[i] = k
        if k == 0:
            disconnected[i] = 0.0
            continue

        failed_idx = rng.choice(m, size=k, replace=False)
        failed_edges = [real_edges[j] for j in failed_idx]
        working.remove_edges_from(failed_edges)
        source_comp = nx.node_connected_component(working, source)
        connected_weight = sum(weights[node] for node in source_comp)
        disconnected[i] = 1.0 - connected_weight / total_weight
        working.add_edges_from((u, v, edge_data[(u, v)]) for u, v in failed_edges)

        if show_progress_every and (i + 1) % show_progress_every == 0:
            elapsed = time.perf_counter() - t0
            print(f"  simulated {i + 1:,}/{n_samples:,} samples in {elapsed:.1f}s")

    mean = float(disconnected.mean())
    std = float(disconnected.std(ddof=1)) if n_samples > 1 else 0.0
    se = std / math.sqrt(n_samples)
    return {
        "method": "stationary_independent_edge_state_monte_carlo",
        "p": p,
        "n_samples": n_samples,
        "seed": seed,
        "mean_disconnected_load_fraction_saidi": mean,
        "std_disconnected_load_fraction": std,
        "standard_error": se,
        "ci95_half_width": 1.96 * se,
        "total_weight_kw": total_weight,
        "real_edges_sampled": m,
        "mean_failed_edges_per_sample": float(failed_edge_counts.mean()),
        "expected_failed_edges_per_sample": m * p,
        "samples_with_no_failed_edges": int((failed_edge_counts == 0).sum()),
        "runtime_seconds": time.perf_counter() - t0,
    }


def write_outputs(metrics: dict) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    json_path = OUT_DIR / "p2u_full_mv_reliability_metrics.json"
    md_path = OUT_DIR / "p2u_full_mv_reliability_metrics.md"
    json_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    topo = metrics["topology"]
    backbone = metrics["backbone"]
    consumers = metrics["consumers"]
    sim = metrics["simulation"]
    lv = metrics["lv_assignment"]
    md_path.write_text(
        "\n".join(
            [
                "# P2U Full MV Reliability Metrics",
                "",
                "## Topology",
                "",
                f"- Nodes: `{topo['nodes']}`",
                f"- Unique MV line edges: `{topo['lines_unique_edges']}`",
                f"- Raw MV line features: `{topo['raw_mv_line_features']}`",
                f"- Normally open tie edges R: `{topo['open_tie_edges_R']}`",
                f"- Physical cycle-rank total ties: `{topo['total_ties_cycle_rank_physical_graph']}`",
                f"- MV switching-device features: `{topo['mv_switching_device_features']}`",
                f"- Closed operating components: `{topo['closed_operating_components']}`",
                f"- Physical components with ties: `{topo['physical_components']}`",
                "",
                "## 2-Edge-Connected Backbone",
                "",
                f"- Source 2-edge-connected backbone nodes: `{backbone['source_2edge_backbone_nodes']}`",
                f"- Nodes in any nontrivial 2-edge-connected component: `{backbone['nodes_in_any_2edge_component']}`",
                f"- Largest 2-edge-connected component nodes: `{backbone['largest_2edge_component_nodes']}`",
                f"- Bridge edges after source contraction: `{backbone['bridge_edges_after_source_contraction']}`",
                "",
                "## Consumers On Source 2-Edge-Connected Backbone",
                "",
                f"- Transformer load nodes: `{consumers['transformer_load_nodes']}`",
                f"- Transformer load nodes 2-edge-connected to source: `{consumers['transformer_load_nodes_2edge_to_source']}`",
                f"- Percent transformer load nodes 2-edge-connected to source: `{consumers['percent_transformer_load_nodes_2edge_to_source']:.3f}%`",
                f"- Demand 2-edge-connected to source: `{consumers['demand_kw_2edge_to_source']:.3f}` kW / `{consumers['demand_kw_total']:.3f}` kW",
                f"- Percent demand 2-edge-connected to source: `{consumers['percent_demand_kw_2edge_to_source']:.3f}%`",
                f"- Customer count 2-edge-connected to source: `{consumers['num_customers_2edge_to_source']:.3f}` / `{consumers['num_customers_total']:.3f}`",
                f"- Percent customers 2-edge-connected to source: `{consumers['percent_customers_2edge_to_source']:.3f}%`",
                "",
                "## LV Assignment",
                "",
                f"- LV consumer points assigned to MV transformers: `{lv['consumer_points_assigned_to_mv']}` / `{lv['consumer_points_total']}`",
                f"- Consumer points assigned by nearest-transformer fallback: `{lv['consumer_points_assigned_by_nearest_transformer_fallback']}`",
                f"- Nearest-transformer fallback mean/max distance: `{lv['nearest_transformer_fallback_mean_distance_m']:.2f}` m / `{lv['nearest_transformer_fallback_max_distance_m']:.2f}` m",
                f"- Assigned demand: `{lv['demand_kw_assigned_to_mv']:.3f}` kW / `{lv['demand_kw_total_raw']:.3f}` kW",
                f"- Transformer nodes with assigned consumers: `{lv['transformer_nodes_with_assigned_consumers']}`",
                "",
                "## SAIDI Monte Carlo",
                "",
                f"- Method: `{sim['method']}`",
                f"- Edge failure probability p: `{sim['p']}`",
                f"- Samples: `{sim['n_samples']}`",
                f"- Seed: `{sim['seed']}`",
                f"- Mean disconnected load fraction SAIDI: `{sim['mean_disconnected_load_fraction_saidi']:.8g}`",
                f"- 95% CI half-width: `{sim['ci95_half_width']:.3g}`",
                f"- Expected failed edges per sample: `{sim['expected_failed_edges_per_sample']:.3f}`",
                f"- Mean failed edges per sample: `{sim['mean_failed_edges_per_sample']:.3f}`",
                f"- Runtime seconds: `{sim['runtime_seconds']:.2f}`",
                "",
                "## Assumptions",
                "",
                "- MV network uses all `NomV` between 1 kV and 40 kV.",
                "- Normally open ties are `Status = 0` MV lines.",
                "- SAIDI is a stationary independent-edge Monte Carlo estimate, not the slower event-loop simulator.",
                "- LV consumer demand is aggregated to MV transformer nodes by LV connected component.",
                "- LV consumers not found in an LV line component are assigned to the nearest distribution transformer.",
                "- 2-connected means 2-edge-connected relative to the contracted source set, because the failure model is edge/line failure.",
            ]
        ),
        encoding="utf-8",
    )
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze full P2U MV physical-network reliability.")
    parser.add_argument("--p", type=float, default=5e-4, help="Independent MV edge failure probability.")
    parser.add_argument("--n-samples", type=int, default=50_000, help="Monte Carlo state samples.")
    parser.add_argument("--seed", type=int, default=20260810, help="Random seed.")
    parser.add_argument("--progress-every", type=int, default=5_000, help="Print progress every N samples; 0 disables.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    graph, mv_rows, switches = build_full_mv_graph()
    lv_assignment = aggregate_lv_consumers_to_transformers(graph)
    sources = graph.graph["sources"]
    if not sources:
        raise ValueError("No P2U MV source nodes were found.")

    source_backbone, backbone = two_edge_backbone_metrics(graph, sources)
    metrics = {
        "dataset": "SMART-DS SFO P2U full MV network",
        "topology": topology_metrics(graph, mv_rows, switches),
        "sources": sources,
        "lv_assignment": lv_assignment,
        "backbone": backbone,
        "consumers": consumer_backbone_metrics(graph, source_backbone),
        "simulation": stationary_mc_saidi(
            graph,
            sources,
            p=args.p,
            n_samples=args.n_samples,
            seed=args.seed,
            show_progress_every=args.progress_every,
        ),
    }
    write_outputs(metrics)


if __name__ == "__main__":
    main()
