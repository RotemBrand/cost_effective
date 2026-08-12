from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.p2u_final_network import (
    FINAL_NETWORK_GPKG,
    FINAL_NETWORK_METADATA,
    OUTPUT_DIR,
    SUPER_SOURCE,
    cycle_rank,
    graph_from_final_tables,
    percent,
    read_final_network_tables,
    source_contracted_graph,
)
from indexes import edge_probs_by_length


SUMMARY_JSON = OUTPUT_DIR / "p2u_final_network_analysis.json"
SUMMARY_MD = OUTPUT_DIR / "p2u_final_network_analysis.md"


def _stationary_terminal_mc(
    graph: nx.Graph,
    sources: set[str],
    *,
    n_samples: int,
    seed: int,
) -> dict:
    rng = np.random.default_rng(seed)
    working = source_contracted_graph(graph, sources)
    if SUPER_SOURCE not in working:
        raise ValueError("source-contracted graph has no source node")

    edges = list(working.edges())
    probs = np.array([float(working.edges[e].get("prob", 0.0)) for e in edges], dtype=float)
    edge_data = {(u, v): working.edges[u, v].copy() for u, v in edges}
    node_weight = {node: float(data.get("weight", 0.0)) for node, data in working.nodes(data=True)}
    edge_weight = {edge: float(working.edges[edge].get("edge_weight", 0.0)) for edge in edges}
    total_weight = float(sum(node_weight.values()) + sum(edge_weight.values()))
    if total_weight <= 0:
        raise ValueError("total demand weight is zero")

    disconnected = np.empty(n_samples, dtype=float)
    failed_counts = np.empty(n_samples, dtype=np.int32)
    t0 = time.perf_counter()
    for i in range(n_samples):
        failed_idx = np.flatnonzero(rng.random(len(edges)) < probs)
        failed_counts[i] = int(len(failed_idx))
        if len(failed_idx):
            failed_edges = [edges[j] for j in failed_idx]
            working.remove_edges_from(failed_edges)
        source_comp = nx.node_connected_component(working, SUPER_SOURCE)
        connected_node_weight = sum(node_weight[node] for node in source_comp)
        connected_edge_weight = sum(
            edge_weight[(u, v)]
            for u, v in working.edges()
            if u in source_comp and v in source_comp
        )
        disconnected[i] = 1.0 - (connected_node_weight + connected_edge_weight) / total_weight
        if len(failed_idx):
            working.add_edges_from((u, v, edge_data[(u, v)]) for u, v in failed_edges)

    mean = float(disconnected.mean())
    std = float(disconnected.std(ddof=1)) if n_samples > 1 else 0.0
    se = std / math.sqrt(n_samples)
    return {
        "method": "terminal_contracted_stationary_monte_carlo",
        "n_samples": int(n_samples),
        "seed": int(seed),
        "mean_disconnected_demand_fraction": mean,
        "std_disconnected_demand_fraction": std,
        "standard_error": se,
        "ci95_half_width": 1.96 * se,
        "samples_with_no_failed_edges": int((failed_counts == 0).sum()),
        "mean_failed_edges_per_sample": float(failed_counts.mean()),
        "runtime_seconds": float(time.perf_counter() - t0),
        "total_weight_kw": total_weight,
        "edge_weight_kw": float(sum(edge_weight.values())),
        "node_weight_kw": float(sum(node_weight.values())),
    }


def analyze_final_network(
    *,
    gpkg_path: Path = FINAL_NETWORK_GPKG,
    metadata_json: Path = FINAL_NETWORK_METADATA,
    p_mean: float = 5e-4,
    n_samples: int = 0,
    seed: int = 20260812,
    output_json: Path = SUMMARY_JSON,
    output_md: Path = SUMMARY_MD,
) -> dict:
    metadata = json.loads(metadata_json.read_text(encoding="utf-8")) if metadata_json.exists() else {}
    if "output_gpkg" in metadata:
        gpkg_path = Path(metadata["output_gpkg"])

    nodes, edges, backbone_edges, tree_edges = read_final_network_tables(gpkg_path)
    graph = graph_from_final_tables(nodes, edges)
    sources = set(nodes[nodes["source_count"].astype(float) > 0]["terminal_id"].astype(str))
    source_graph = source_contracted_graph(graph, sources)

    edge_probs, failure_rate = edge_probs_by_length(graph, p=p_mean, mode="mean", length_attr="length_m")
    for (u, v), prob in edge_probs.items():
        graph.edges[u, v]["prob"] = prob

    source_graph_bridges = list(nx.bridges(source_graph)) if nx.is_connected(source_graph) else []
    node_transformers = int(nodes["transformer_count"].sum())
    edge_transformers = int(edges["edge_transformer_count"].sum())
    represented_transformers = node_transformers + edge_transformers
    total_transformers = int(metadata.get("total_transformers_in_data", represented_transformers))

    total_size_kva = float(nodes["size_kva"].sum() + edges["edge_size_kva"].sum())
    backbone_size_kva = float(
        nodes[nodes["terminal_class"].isin(["backbone_boundary", "source"])]["size_kva"].sum()
        + backbone_edges["edge_size_kva"].sum()
    )
    total_demand_kw = float(nodes["demand_kw"].sum() + edges["edge_demand_kw"].sum())
    backbone_demand_kw = float(
        nodes[nodes["terminal_class"].isin(["backbone_boundary", "source"])]["demand_kw"].sum()
        + backbone_edges["edge_demand_kw"].sum()
    )
    tree_demand_kw = float(nodes[nodes["terminal_class"] == "tree_terminal"]["demand_kw"].sum())
    backbone_transformers = int(
        nodes[nodes["terminal_class"].isin(["backbone_boundary", "source"])]["transformer_count"].sum()
        + backbone_edges["edge_transformer_count"].sum()
    )
    tree_transformers = int(nodes[nodes["terminal_class"] == "tree_terminal"]["transformer_count"].sum())

    simulation = None
    if n_samples > 0:
        simulation = _stationary_terminal_mc(graph, sources, n_samples=n_samples, seed=seed)

    summary = {
        "input_gpkg": str(gpkg_path),
        "metadata_json": str(metadata_json),
        "final_graph_definition": metadata.get("final_graph_definition", ""),
        "nodes": int(graph.number_of_nodes()),
        "edges": int(graph.number_of_edges()),
        "connected_components": int(nx.number_connected_components(graph)),
        "source_nodes": int(len(sources)),
        "source_contracted_nodes": int(source_graph.number_of_nodes()),
        "source_contracted_edges": int(source_graph.number_of_edges()),
        "source_contracted_components": int(nx.number_connected_components(source_graph)),
        "source_contracted_cycle_rank_R": int(cycle_rank(source_graph)),
        "source_contracted_bridge_edges": int(len(source_graph_bridges)),
        "source_contracted_is_connected": bool(nx.is_connected(source_graph)),
        "backbone_edges": int(len(backbone_edges)),
        "tree_attachment_edges": int(len(tree_edges)),
        "backbone_boundary_transformer_nodes": int((nodes["terminal_class"] == "backbone_boundary").sum()),
        "tree_transformer_terminal_nodes": int((nodes["terminal_class"] == "tree_terminal").sum()),
        "total_transformers_in_data": total_transformers,
        "represented_transformers": represented_transformers,
        "node_transformers": node_transformers,
        "edge_load_transformers": edge_transformers,
        "backbone_transformers_including_edge_load": backbone_transformers,
        "tree_transformers": tree_transformers,
        "percent_transformers_on_backbone": percent(backbone_transformers, represented_transformers),
        "total_size_kva": total_size_kva,
        "backbone_size_kva_including_edge_load": backbone_size_kva,
        "percent_size_kva_on_backbone": percent(backbone_size_kva, total_size_kva),
        "total_demand_kw": total_demand_kw,
        "backbone_demand_kw_including_edge_load": backbone_demand_kw,
        "tree_demand_kw": tree_demand_kw,
        "percent_demand_on_backbone": percent(backbone_demand_kw, total_demand_kw),
        "backbone_length_m": float(backbone_edges["length_m"].sum()),
        "tree_attachment_terminal_length_sum_m": float(tree_edges["length_m"].sum()) if len(tree_edges) else 0.0,
        "tree_attachment_physical_road_union_length_m": float(
            metadata.get("tree_attachment_physical_road_union_length_m", 0.0)
        ),
        "total_physical_length_m": float(backbone_edges["length_m"].sum())
        + float(metadata.get("tree_attachment_physical_road_union_length_m", 0.0)),
        "mean_tree_terminal_distance_to_backbone_m": float(tree_edges["length_m"].mean()) if len(tree_edges) else 0.0,
        "max_tree_terminal_distance_to_backbone_m": float(tree_edges["length_m"].max()) if len(tree_edges) else 0.0,
        "p_mean_target": float(p_mean),
        "length_failure_rate_per_m": float(failure_rate),
        "simulation": simulation,
    }
    output_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    output_md.write_text(summary_markdown(summary), encoding="utf-8")
    return summary


def summary_markdown(summary: dict) -> str:
    sim = summary.get("simulation")
    sim_lines = ["- Monte Carlo: `not run`"]
    if sim is not None:
        sim_lines = [
            f"- Monte Carlo samples: `{sim['n_samples']}`",
            f"- Mean disconnected demand fraction: `{sim['mean_disconnected_demand_fraction']:.8g}`",
            f"- 95% CI half-width: `{sim['ci95_half_width']:.3g}`",
            f"- Mean failed edges/sample: `{sim['mean_failed_edges_per_sample']:.3f}`",
            f"- Runtime: `{sim['runtime_seconds']:.2f}` s",
            "- Interpretation: switch-section graph; edge loads represent contracted transformer chains, and zero-load street branch nodes may represent shared radial road branches.",
        ]
    return "\n".join(
        [
            "# P2U Final Optimized Network Analysis",
            "",
            "The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.",
            "",
            "## Topology",
            "",
            f"- Final graph nodes/edges: `{summary['nodes']}` / `{summary['edges']}`",
            f"- Source-contracted nodes/edges: `{summary['source_contracted_nodes']}` / `{summary['source_contracted_edges']}`",
            f"- Source-contracted cycle rank R: `{summary['source_contracted_cycle_rank_R']}`",
            f"- Source-contracted bridge edges: `{summary['source_contracted_bridge_edges']}`",
            f"- Backbone edges: `{summary['backbone_edges']}`",
            f"- Tree attachment edges: `{summary['tree_attachment_edges']}`",
            f"- Sources: `{summary['source_nodes']}`",
            "",
            "## Transformers And Load",
            "",
            f"- Represented transformers: `{summary['represented_transformers']}` / `{summary['total_transformers_in_data']}`",
            f"- Node transformers: `{summary['node_transformers']}`",
            f"- Contracted edge-load transformers: `{summary['edge_load_transformers']}`",
            f"- Transformers on backbone, including edge load: `{summary['backbone_transformers_including_edge_load']}`",
            f"- Tree transformers: `{summary['tree_transformers']}`",
            f"- Percent transformers on backbone: `{summary['percent_transformers_on_backbone']:.3f}%`",
            f"- Capacity on backbone: `{summary['backbone_size_kva_including_edge_load']:.1f}` kVA / `{summary['total_size_kva']:.1f}` kVA",
            f"- Demand on backbone: `{summary['backbone_demand_kw_including_edge_load']:.3f}` kW / `{summary['total_demand_kw']:.3f}` kW",
            f"- Percent demand on backbone: `{summary['percent_demand_on_backbone']:.3f}%`",
            "",
            "## Length",
            "",
            f"- Backbone length: `{summary['backbone_length_m']:.3f}` m",
            f"- Tree physical road-union length: `{summary['tree_attachment_physical_road_union_length_m']:.3f}` m",
            f"- Total physical length: `{summary['total_physical_length_m']:.3f}` m",
            f"- Terminal attachment distance sum: `{summary['tree_attachment_terminal_length_sum_m']:.3f}` m",
            f"- Mean/max terminal distance to backbone: `{summary['mean_tree_terminal_distance_to_backbone_m']:.3f}` m / `{summary['max_tree_terminal_distance_to_backbone_m']:.3f}` m",
            "",
            "## Failure Probabilities",
            "",
            f"- Target mean edge probability: `{summary['p_mean_target']}`",
            f"- Fitted failure rate per meter: `{summary['length_failure_rate_per_m']:.8g}`",
            "",
            "## Reliability",
            "",
            *sim_lines,
            "",
            "## Notes",
            "",
            "- The tree attachment edge length is the road shortest-path distance to the selected backbone.",
            "- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.",
        ]
    ) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze the final P2U optimized network.")
    parser.add_argument("--gpkg", type=Path, default=FINAL_NETWORK_GPKG)
    parser.add_argument("--metadata-json", type=Path, default=FINAL_NETWORK_METADATA)
    parser.add_argument("--p-mean", type=float, default=5e-4)
    parser.add_argument("--n-samples", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--output-json", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--output-md", type=Path, default=SUMMARY_MD)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = analyze_final_network(
        gpkg_path=args.gpkg,
        metadata_json=args.metadata_json,
        p_mean=args.p_mean,
        n_samples=args.n_samples,
        seed=args.seed,
        output_json=args.output_json,
        output_md=args.output_md,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
