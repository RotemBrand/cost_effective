from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from optimal_network import optimize_rel_weight_ratio  # noqa: E402
from figures.optimal_sfo.run_p2u_euclidean_equal_chains import (  # noqa: E402
    OUTPUT_DIR,
    OUTPUT_GPKG as DEFAULT_INPUT_GPKG,
    cycle_rank,
    graph_to_geodataframes,
    internal_chain_risk_proxy,
    refresh_edge_attributes,
    total_graph_length,
    write_outputs,
)


DEFAULT_OUTPUT_GPKG = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_localfix_3857.gpkg"
DEFAULT_OUTPUT_QGS = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_localfix.qgs"
DEFAULT_SUMMARY_JSON = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_localfix_summary.json"
DEFAULT_SUMMARY_MD = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_localfix_summary.md"


def read_euclidean_graph(gpkg: Path) -> tuple[nx.Graph, gpd.GeoDataFrame]:
    nodes = gpd.read_file(gpkg, layer="euclidean_nodes")
    edges = gpd.read_file(gpkg, layer="euclidean_edges")
    graph = nx.Graph()
    for _, row in nodes.iterrows():
        node = int(row.point_index)
        graph.add_node(
            node,
            terminal_id=str(row.terminal_id),
            kind=str(row.kind),
            source_count=int(row.source_count),
            transformer_count=int(row.transformer_count),
            size_kva=float(row.size_kva),
            nominal_voltage_kv=float(row.nominal_voltage_kv),
            pos=np.asarray([float(row.geometry.x), float(row.geometry.y)], dtype=float),
        )
    for _, row in edges.iterrows():
        graph.add_edge(int(row.u_index), int(row.v_index))
    refresh_edge_attributes(graph)
    return graph, nodes


def apply_local_fixes(
    *,
    input_gpkg: Path,
    max_changes: int,
    max_risk_gain: float,
    output_gpkg: Path,
    output_qgs: Path,
    summary_json: Path,
    summary_md: Path,
    debug: bool,
) -> dict:
    graph, terminals = read_euclidean_graph(input_gpkg)
    before_length = total_graph_length(graph)
    before_risk = internal_chain_risk_proxy(graph)
    before_r = cycle_rank(graph)
    before_bridges = len(list(nx.bridges(graph))) if nx.is_connected(graph) else None

    t0 = time.perf_counter()
    fixed = optimize_rel_weight_ratio(
        graph,
        max_risk_gain=max_risk_gain,
        max_changes=max_changes,
        source=None,
        debug=debug,
    )
    refresh_edge_attributes(fixed)
    runtime = time.perf_counter() - t0

    after_r = cycle_rank(fixed)
    after_bridges = len(list(nx.bridges(fixed))) if nx.is_connected(fixed) else None
    if after_r != before_r:
        raise RuntimeError(f"local fixes changed cycle rank from {before_r} to {after_r}")
    if before_bridges == 0 and after_bridges != 0:
        raise RuntimeError(f"local fixes introduced {after_bridges} bridges")

    fixed.graph["local_fix_summary"] = {
        "enabled": True,
        "input_gpkg": str(input_gpkg),
        "max_changes": int(max_changes),
        "max_risk_gain": float(max_risk_gain),
        "runtime_s": float(runtime),
        "length_before_m": float(before_length),
        "length_after_m": float(total_graph_length(fixed)),
        "length_delta_m": float(total_graph_length(fixed) - before_length),
        "length_relative_change": float(total_graph_length(fixed) / before_length - 1.0) if before_length else 0.0,
        "internal_chain_risk_before": float(before_risk),
        "internal_chain_risk_after": float(internal_chain_risk_proxy(fixed)),
        "cycle_rank_before": int(before_r),
        "cycle_rank_after": int(after_r),
        "bridges_before": before_bridges,
        "bridges_after": after_bridges,
    }

    tables = graph_to_geodataframes(fixed, terminals)
    args = argparse.Namespace(
        input_gpkg=input_gpkg,
        seed=0,
        kmeans_max_iter=0,
        strc_n_init_iters=0,
        strc_exact_vertices=False,
        strc_trip_nearest_vertices=None,
        chain_n_init_iters=0,
    )
    summary = write_outputs(
        tables=tables,
        graph=fixed,
        redundancy=before_r,
        runtime_s=runtime,
        output_gpkg=output_gpkg,
        output_qgs=output_qgs,
        summary_json=summary_json,
        summary_md=summary_md,
        args=args,
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Apply greedy local cost fixes to an existing P2U Euclidean network.")
    parser.add_argument("--input-gpkg", type=Path, default=DEFAULT_INPUT_GPKG)
    parser.add_argument("--max-changes", type=int, default=50)
    parser.add_argument("--max-risk-gain", type=float, default=0.05)
    parser.add_argument("--output-gpkg", type=Path, default=DEFAULT_OUTPUT_GPKG)
    parser.add_argument("--output-qgs", type=Path, default=DEFAULT_OUTPUT_QGS)
    parser.add_argument("--summary-json", type=Path, default=DEFAULT_SUMMARY_JSON)
    parser.add_argument("--summary-md", type=Path, default=DEFAULT_SUMMARY_MD)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    summary = apply_local_fixes(
        input_gpkg=args.input_gpkg,
        max_changes=args.max_changes,
        max_risk_gain=args.max_risk_gain,
        output_gpkg=args.output_gpkg,
        output_qgs=args.output_qgs,
        summary_json=args.summary_json,
        summary_md=args.summary_md,
        debug=args.debug,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
