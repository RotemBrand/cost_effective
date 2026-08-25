from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.analyze_p2u_final_network import analyze_final_network
from figures.optimal_sfo.build_p2u_final_network import build_and_write_final_network
from figures.optimal_sfo.create_p2u_final_network_qgis import DEFAULT_QGS, write_qgis_project
from figures.optimal_sfo.prepare_p2u_corridor_network import build_corridor_outputs
from figures.optimal_sfo.run_p2u_ilp_2edge import load_ilp_graph, solve_min_2edge, write_solution


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the modular P2U optimal-network pipeline.")
    parser.add_argument("--skip-preparation", action="store_true")
    parser.add_argument("--skip-ilp", action="store_true")
    parser.add_argument("--skip-final-network", action="store_true")
    parser.add_argument("--skip-qgis", action="store_true")
    parser.add_argument("--skip-analysis", action="store_true")
    parser.add_argument("--ilp-time-limit", type=float, default=600.0)
    parser.add_argument("--ilp-mip-gap", type=float, default=0.02)
    parser.add_argument("--ilp-threads", type=int, default=0)
    parser.add_argument("--ilp-max-cut-rounds", type=int, default=100)
    parser.add_argument("--ilp-max-redundancy", type=int, default=None)
    parser.add_argument("--ilp-cut-mode", choices=["iterative", "callback"], default="callback")
    parser.add_argument(
        "--ilp-objective-mode",
        choices=["min_length", "max_chain_demand_then_min_length"],
        default="max_chain_demand_then_min_length",
    )
    parser.add_argument("--ilp-coverage-attr", default="edge_size_kva")
    parser.add_argument("--ilp-coverage-tolerance", type=float, default=1e-6)
    parser.add_argument("--tree-mode", choices=["street_forest", "star"], default="street_forest")
    parser.add_argument("--p-mean", type=float, default=5e-4)
    parser.add_argument("--n-samples", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20260812)
    return parser.parse_args()


def run_ilp_stage(args: argparse.Namespace) -> dict:
    graph, edges, transformer_nodes, source_nodes, original_source_incident = load_ilp_graph()
    solution, summary = solve_min_2edge(
        graph,
        original_source_incident=original_source_incident,
        redundancy=None,
        max_redundancy=args.ilp_max_redundancy,
        time_limit=args.ilp_time_limit,
        mip_gap=args.ilp_mip_gap,
        threads=args.ilp_threads,
        max_cut_rounds=args.ilp_max_cut_rounds,
        cut_mode=args.ilp_cut_mode,
        objective_mode=args.ilp_objective_mode,
        coverage_attr=args.ilp_coverage_attr,
        coverage_tolerance=args.ilp_coverage_tolerance,
    )
    write_solution(solution, summary, edges, transformer_nodes, source_nodes)
    return summary


def main() -> None:
    args = parse_args()
    results = {}

    if not args.skip_preparation:
        results["preparation"] = build_corridor_outputs()

    if not args.skip_ilp:
        results["ilp"] = run_ilp_stage(args)

    final_metadata = None
    if not args.skip_final_network:
        final_metadata = build_and_write_final_network(tree_mode=args.tree_mode)
        results["final_network"] = final_metadata

    if not args.skip_qgis:
        gpkg = None if final_metadata is None else final_metadata.get("output_gpkg")
        results["qgis"] = write_qgis_project(gpkg_path=gpkg, qgs_path=DEFAULT_QGS) if gpkg else write_qgis_project()

    if not args.skip_analysis:
        results["analysis"] = analyze_final_network(
            p_mean=args.p_mean,
            n_samples=args.n_samples,
            seed=args.seed,
        )

    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
