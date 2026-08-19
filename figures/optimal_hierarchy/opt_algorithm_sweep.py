from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import networkx as nx
import pandas as pd
from shapely.geometry import LineString

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_hierarchy.redundancy_sweep import (  # noqa: E402
    SweepParameters,
    _case_row,
    compute_original_reference,
)
from figures.optimal_sfo.analyze_p2u_final_network import analyze_final_network  # noqa: E402
from figures.optimal_sfo.compare_p2u_old_new_reliability import (  # noqa: E402
    _analyze_one,
    _prepare_new_graph,
)
from figures.optimal_sfo.p2u_final_network import build_and_write_final_network  # noqa: E402
from figures.optimal_sfo.run_p2u_euclidean_equal_chains import (  # noqa: E402
    cycle_rank,
    run_equal_chain_algorithm,
)


CORRIDOR_GPKG = ROOT / "outputs" / "optimal_sfo" / "p2u_terminal_corridors_road2_k10_3857.gpkg"
OUTPUT_DIR = ROOT / "outputs" / "optimal_hierarchy" / "opt"


@dataclass(frozen=True)
class OptParameters:
    p_mean: float = 5e-4
    seed: int = 20260812
    kmeans_max_iter: int = 3
    strc_n_init_iters: int = 1
    strc_exact_vertices: bool = True
    strc_trip_nearest_vertices: int | None = 24
    chain_n_init_iters: int = 1
    local_fix_max_changes: int = 0
    local_fix_max_risk_gain: float = 0.2
    generalized_method: str = "projection"
    tree_mode: str = "street_forest"
    target_source_contracted_r: bool = True


def try_dir(try_name: str) -> Path:
    if try_name in {"", ".", "opt"}:
        return OUTPUT_DIR
    return OUTPUT_DIR / try_name


def case_dir(try_name: str, r_value: int) -> Path:
    return try_dir(try_name) / f"R{r_value:03d}"


def case_paths(try_name: str, r_value: int) -> dict[str, Path]:
    directory = case_dir(try_name, r_value)
    token = f"R{r_value}"
    return {
        "dir": directory,
        "backbone_gpkg": directory / f"p2u_opt_backbone_{token}_3857.gpkg",
        "backbone_json": directory / f"p2u_opt_backbone_{token}_summary.json",
        "backbone_md": directory / f"p2u_opt_backbone_{token}_summary.md",
        "final_gpkg": directory / f"p2u_opt_final_network_{token}_streetforest_3857.gpkg",
        "final_metadata": directory / f"p2u_opt_final_network_{token}_streetforest_metadata.json",
        "topology_json": directory / f"p2u_opt_final_network_{token}_topology.json",
        "topology_md": directory / f"p2u_opt_final_network_{token}_topology.md",
        "risk_json": directory / f"p2u_opt_final_network_{token}_risk.json",
    }


def load_backbone_terminals(corridor_gpkg: Path) -> gpd.GeoDataFrame:
    transformers = gpd.read_file(corridor_gpkg, layer="ilp_transformer_nodes")
    sources = gpd.read_file(corridor_gpkg, layer="ilp_source_nodes")
    if transformers.crs != sources.crs:
        sources = sources.to_crs(transformers.crs)
    terminals = gpd.GeoDataFrame(
        pd.concat([transformers, sources], ignore_index=True),
        geometry="geometry",
        crs=transformers.crs,
    )
    terminals = terminals.drop_duplicates(subset=["terminal_id"], keep="first").reset_index(drop=True)
    terminals["point_index"] = terminals.index.astype(int)
    if "nominal_voltage_kv" not in terminals.columns:
        terminals["nominal_voltage_kv"] = 12.47
    return terminals


def build_opt_backbone(
    *,
    terminals: gpd.GeoDataFrame,
    physical_redundancy: int,
    params: OptParameters,
    checkpoint_dir: Path,
    debug: bool,
) -> nx.Graph:
    graph = run_equal_chain_algorithm(
        terminals,
        redundancy=physical_redundancy,
        seed=params.seed,
        kmeans_max_iter=params.kmeans_max_iter,
        strc_n_init_iters=params.strc_n_init_iters,
        strc_exact_vertices=params.strc_exact_vertices,
        strc_trip_nearest_vertices=params.strc_trip_nearest_vertices,
        chain_n_init_iters=params.chain_n_init_iters,
        local_fix_max_changes=params.local_fix_max_changes,
        local_fix_max_risk_gain=params.local_fix_max_risk_gain,
        debug=debug,
        checkpoint_dir=checkpoint_dir,
    )
    if not nx.is_connected(graph):
        raise RuntimeError("OPT backbone graph is not connected")
    bridges = list(nx.bridges(graph))
    if bridges:
        raise RuntimeError(f"OPT backbone graph has {len(bridges)} bridges")
    achieved_r = cycle_rank(graph)
    if achieved_r != physical_redundancy:
        raise RuntimeError(f"OPT backbone achieved physical R={achieved_r}, expected R={physical_redundancy}")
    return graph


def write_opt_backbone_solution(
    *,
    graph: nx.Graph,
    terminals: gpd.GeoDataFrame,
    requested_source_contracted_r: int,
    physical_redundancy: int,
    params: OptParameters,
    runtime_s: float,
    paths: dict[str, Path],
) -> dict[str, Any]:
    paths["dir"].mkdir(parents=True, exist_ok=True)
    terminal_rows = terminals.set_index("point_index")

    edge_rows = []
    for edge_id, (u, v, data) in enumerate(graph.edges(data=True)):
        row_u = terminal_rows.loc[int(u)]
        row_v = terminal_rows.loc[int(v)]
        terminal_a = str(row_u.terminal_id)
        terminal_b = str(row_v.terminal_id)
        road_node_a = str(row_u.road_node)
        road_node_b = str(row_v.road_node)
        pos_u = graph.nodes[u]["pos"]
        pos_v = graph.nodes[v]["pos"]
        edge_rows.append(
            {
                "edge_id": int(edge_id),
                "terminal_a": terminal_a,
                "terminal_b": terminal_b,
                "road_node_a": road_node_a,
                "road_node_b": road_node_b,
                "terminal_kind_a": str(row_u.kind),
                "terminal_kind_b": str(row_v.kind),
                "edge_transformer_count": 0,
                "edge_source_count": 0,
                "edge_size_kva": 0.0,
                "chain_terminal_count": 2,
                "internal_terminal_count": 0,
                "length_m": float(data["length_m"]),
                "road_node_count": 0,
                "terminal_sequence": f"{terminal_a}|{terminal_b}",
                "is_tie": bool(data.get("is_tie", False)),
                "is_switch": bool(data.get("is_switch", False)),
                "normally_closed": bool(data.get("normally_closed", True)),
                "geometry": LineString([pos_u, pos_v]),
            }
        )
    edges = gpd.GeoDataFrame(edge_rows, geometry="geometry", crs=terminals.crs)
    selected_ids = {int(node) for node in graph.nodes}
    nodes = terminals[terminals["point_index"].isin(selected_ids)].copy()
    transformer_nodes = nodes[nodes["kind"] == "transformer"].copy()
    source_nodes = nodes[nodes["kind"] == "source"].copy()

    if paths["backbone_gpkg"].exists():
        paths["backbone_gpkg"].unlink()
    edges.to_file(paths["backbone_gpkg"], layer="solution_edges", driver="GPKG")
    source_mask = nodes["source_count"].astype(float) > 0
    source_nodes = nodes[source_mask].copy()
    transformer_nodes = nodes[~source_mask & (nodes["transformer_count"].astype(float) > 0)].copy()
    transformer_nodes.to_file(paths["backbone_gpkg"], layer="solution_transformer_nodes", driver="GPKG")
    source_nodes.to_file(paths["backbone_gpkg"], layer="solution_source_nodes", driver="GPKG")

    lengths = [float(data["length_m"]) for _, _, data in graph.edges(data=True)]
    summary = {
        "algorithm": "OPT_equal_chain_backbone",
        "output_gpkg": str(paths["backbone_gpkg"]),
        "redundancy_constraint": int(requested_source_contracted_r),
        "physical_redundancy_constraint": int(physical_redundancy),
        "target_source_contracted_r": bool(params.target_source_contracted_r),
        "max_redundancy_constraint": None,
        "status_name": "OPT_CONSTRUCTED",
        "stop_reason": "opt_algorithm_completed",
        "runtime_s": float(runtime_s),
        "seed": int(params.seed),
        "kmeans_max_iter": int(params.kmeans_max_iter),
        "strc_n_init_iters": int(params.strc_n_init_iters),
        "strc_exact_vertices": bool(params.strc_exact_vertices),
        "strc_trip_nearest_vertices": params.strc_trip_nearest_vertices,
        "chain_n_init_iters": int(params.chain_n_init_iters),
        "local_fix": graph.graph.get("local_fix_summary", {"enabled": False}),
        "input_nodes": int(len(terminals)),
        "solution_nodes": int(graph.number_of_nodes()),
        "solution_edges": int(graph.number_of_edges()),
        "solution_cycle_rank": int(cycle_rank(graph)),
        "solution_connected": bool(nx.is_connected(graph)),
        "solution_bridge_count": int(len(list(nx.bridges(graph)))),
        "solution_is_2edge_connected": bool(nx.is_connected(graph) and len(list(nx.bridges(graph))) == 0),
        "solution_transformer_nodes": int(len(transformer_nodes)),
        "solution_source_nodes": int(len(source_nodes)),
        "objective_length_m": float(sum(lengths)),
        "mean_edge_length_m": float(sum(lengths) / len(lengths)) if lengths else 0.0,
        "max_edge_length_m": float(max(lengths)) if lengths else 0.0,
        "tie_edges": int(sum(1 for _, _, data in graph.edges(data=True) if data.get("is_tie", False))),
        "tree_mode": params.tree_mode,
    }
    paths["backbone_json"].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    paths["backbone_md"].write_text(_backbone_markdown(summary), encoding="utf-8")
    return summary


def _backbone_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# P2U OPT Backbone",
        "",
        "This backbone replaces the road-corridor ILP with the OPT equal-chain construction on the same selected backbone terminal set.",
        "",
        f"- Status: `{summary['status_name']}`",
        f"- Requested R: `{summary['redundancy_constraint']}`",
        f"- Physical OPT R: `{summary['physical_redundancy_constraint']}`",
        f"- Achieved physical R: `{summary['solution_cycle_rank']}`",
        f"- Nodes/edges: `{summary['solution_nodes']}` / `{summary['solution_edges']}`",
        f"- Sources/transformers: `{summary['solution_source_nodes']}` / `{summary['solution_transformer_nodes']}`",
        f"- 2-edge-connected: `{summary['solution_is_2edge_connected']}`",
        f"- Tie edges: `{summary['tie_edges']}`",
        f"- Euclidean backbone length: `{summary['objective_length_m']:.3f}` m",
        f"- Mean/max edge length: `{summary['mean_edge_length_m']:.3f}` m / `{summary['max_edge_length_m']:.3f}` m",
        f"- Runtime: `{summary['runtime_s']:.2f}` s",
        "",
        "The final hierarchy still attaches non-backbone transformers through the street-forest code.",
    ]
    return "\n".join(lines) + "\n"


def analyze_opt_case(
    *,
    r_value: int,
    params: OptParameters,
    paths: dict[str, Path],
    reuse_existing: bool,
) -> dict[str, Any]:
    if reuse_existing and paths["risk_json"].exists() and paths["topology_json"].exists():
        return json.loads(paths["risk_json"].read_text(encoding="utf-8"))

    topology_t0 = time.perf_counter()
    topology = analyze_final_network(
        gpkg_path=paths["final_gpkg"],
        metadata_json=paths["final_metadata"],
        p_mean=params.p_mean,
        output_json=paths["topology_json"],
        output_md=paths["topology_md"],
    )
    topology["stage_runtime_s"] = time.perf_counter() - topology_t0
    paths["topology_json"].write_text(json.dumps(topology, indent=2), encoding="utf-8")

    graph, sources, extra = _prepare_new_graph(paths["final_gpkg"], paths["final_metadata"])
    extra["topology_summary"] = topology
    risk = _analyze_one(
        name=f"opt_hierarchy_R{r_value}",
        graph=graph,
        sources=sources,
        extra=extra,
        p_mean=params.p_mean,
        generalized=True,
        generalized_method=params.generalized_method,
    )
    risk["topology_summary"] = topology
    paths["risk_json"].write_text(json.dumps(risk, indent=2), encoding="utf-8")
    return risk


def run_case(
    *,
    try_name: str,
    r_value: int,
    terminals: gpd.GeoDataFrame,
    params: OptParameters,
    reuse_existing: bool,
    debug: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    paths = case_paths(try_name, r_value)
    paths["dir"].mkdir(parents=True, exist_ok=True)
    if reuse_existing and paths["backbone_json"].exists() and paths["backbone_gpkg"].exists():
        backbone_summary = json.loads(paths["backbone_json"].read_text(encoding="utf-8"))
    else:
        n_sources = int((terminals["source_count"].astype(float) > 0).sum())
        physical_redundancy = r_value - (n_sources - 1) if params.target_source_contracted_r else r_value
        if physical_redundancy <= 0:
            raise ValueError(
                f"requested source-contracted R={r_value} is too small for {n_sources} sources; "
                f"physical redundancy would be {physical_redundancy}"
            )
        t0 = time.perf_counter()
        graph = build_opt_backbone(
            terminals=terminals,
            physical_redundancy=physical_redundancy,
            params=params,
            checkpoint_dir=paths["dir"],
            debug=debug,
        )
        backbone_summary = write_opt_backbone_solution(
            graph=graph,
            terminals=terminals,
            requested_source_contracted_r=r_value,
            physical_redundancy=physical_redundancy,
            params=params,
            runtime_s=time.perf_counter() - t0,
            paths=paths,
        )

    if not (reuse_existing and paths["final_gpkg"].exists() and paths["final_metadata"].exists()):
        build_and_write_final_network(
            backbone_summary=paths["backbone_json"],
            output_gpkg=paths["final_gpkg"],
            metadata_json=paths["final_metadata"],
            tree_mode=params.tree_mode,
        )
    risk = analyze_opt_case(r_value=r_value, params=params, paths=paths, reuse_existing=reuse_existing)
    return backbone_summary, risk


def write_summary(
    *,
    try_name: str,
    rows: list[dict[str, Any]],
    stage_log: dict[str, Any],
    params: OptParameters,
) -> Path:
    out_dir = try_dir(try_name)
    out_dir.mkdir(parents=True, exist_ok=True)
    table_csv = out_dir / "p2u_opt_hierarchy_table.csv"
    table_json = out_dir / "p2u_opt_hierarchy_table.json"
    manifest = out_dir / "p2u_opt_hierarchy.md"
    pd.DataFrame(rows).to_csv(table_csv, index=False, lineterminator="\n")
    table_json.write_text(json.dumps({"rows": rows, "stage_log": stage_log}, indent=2), encoding="utf-8")

    display_cols = [
        "network",
        "r_request",
        "length_km",
        "r_theory",
        "z_w",
        "z_r",
        "z_f",
        "z_f_p",
        "z_f_p2",
        "risk_total",
        "risk_o_p",
        "risk_o_p2",
        "risk_tree",
        "risk_section",
        "risk_internal",
        "risk_structural",
        "bridges",
        "generalized_chains",
        "gen_lambda_mean_km",
        "gen_lambda_sigma_over_mean",
        "gen_lambda_max_km",
        "ilp_runtime_s",
        "decomposition_runtime_s",
    ]
    lines = [
        f"# P2U OPT Hierarchy: {try_name}",
        "",
        "## Results",
        "",
        "This experiment replaces the road-corridor ILP backbone with the OPT equal-chain construction on the same backbone-terminal set. The final graph still uses the street-forest attachment for non-backbone transformers.",
        "",
        "The OPT backbone is Euclidean, so the backbone cost is not road-constrained. The tree attachment cost remains street-based.",
        "",
        "When `target_source_contracted_r=true`, the OPT physical redundancy is reduced by `n_sources - 1` so the final source-contracted theory graph has the requested `R`.",
        "",
        _markdown_table(pd.DataFrame(rows)[display_cols]),
        "",
        "## Algorithm",
        "",
        "1. Load the same terminal set used by the road-corridor ILP backbone.",
        "2. Run balanced chain clustering and OPT structure construction for the requested `R`.",
        "3. Export the OPT graph as a backbone solution with straight-line Euclidean edges.",
        "4. Reuse the final-network builder to attach all remaining transformers as a street forest.",
        "5. Run the deterministic switch-aware reliability decomposition with generalized chains.",
        "",
        "## Implementation",
        "",
        f"- Script: `figures/optimal_hierarchy/opt_algorithm_sweep.py`",
        f"- Output table: `{table_csv.relative_to(ROOT)}`",
        f"- Output JSON: `{table_json.relative_to(ROOT)}`",
        "",
        "Parameters:",
        "",
        "```json",
        json.dumps(params.__dict__, indent=2),
        "```",
        "",
        "Stage log:",
        "",
        "```json",
        json.dumps(stage_log, indent=2),
        "```",
        "",
    ]
    manifest.write_text("\n".join(lines), encoding="utf-8")
    return manifest


def _format_md_value(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, float):
        if abs(value - round(value)) < 1e-12:
            return str(int(round(value)))
        mantissa, exponent = f"{value:.2e}".split("e")
        return f"{mantissa}e{int(exponent)}"
    return str(value)


def _markdown_table(df: pd.DataFrame) -> str:
    headers = list(df.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(_format_md_value(row[col]) for col in headers) + " |")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run P2U hierarchy cases with an OPT Euclidean backbone.")
    parser.add_argument("--corridor-gpkg", type=Path, default=CORRIDOR_GPKG)
    parser.add_argument(
        "--try-name",
        default="opt",
        help="Output subfolder under outputs/optimal_hierarchy/opt. Use 'opt' for the method folder itself.",
    )
    parser.add_argument("--r-values", type=int, nargs="+", default=[50])
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--kmeans-max-iter", type=int, default=3)
    parser.add_argument("--strc-n-init-iters", type=int, default=1)
    parser.add_argument("--strc-exact-vertices", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--strc-trip-nearest-vertices", type=int, default=24)
    parser.add_argument("--chain-n-init-iters", type=int, default=1)
    parser.add_argument("--local-fix-max-changes", type=int, default=0)
    parser.add_argument("--local-fix-max-risk-gain", type=float, default=0.2)
    parser.add_argument("--p-mean", type=float, default=5e-4)
    parser.add_argument("--generalized-method", choices=["projection", "networkx"], default="projection")
    parser.add_argument("--tree-mode", choices=["street_forest", "star"], default="street_forest")
    parser.add_argument("--target-source-contracted-r", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = OptParameters(
        p_mean=args.p_mean,
        seed=args.seed,
        kmeans_max_iter=args.kmeans_max_iter,
        strc_n_init_iters=args.strc_n_init_iters,
        strc_exact_vertices=args.strc_exact_vertices,
        strc_trip_nearest_vertices=args.strc_trip_nearest_vertices,
        chain_n_init_iters=args.chain_n_init_iters,
        local_fix_max_changes=args.local_fix_max_changes,
        local_fix_max_risk_gain=args.local_fix_max_risk_gain,
        generalized_method=args.generalized_method,
        tree_mode=args.tree_mode,
        target_source_contracted_r=args.target_source_contracted_r,
    )

    t0 = time.perf_counter()
    terminals = load_backbone_terminals(args.corridor_gpkg)
    stage_log: dict[str, Any] = {
        "parameters": params.__dict__,
        "corridor_gpkg": str(args.corridor_gpkg),
        "backbone_terminal_count": int(len(terminals)),
        "cases": {},
    }
    original_params = SweepParameters(p_mean=params.p_mean, generalized_method=params.generalized_method)
    original = compute_original_reference(original_params, reuse_existing=True)
    rows = [
        _case_row(
            r_request=None,
            label="Original P2U MV",
            analysis=original,
            original=original,
            backbone_summary=None,
            is_original=True,
        )
    ]

    for r_value in args.r_values:
        case_t0 = time.perf_counter()
        backbone_summary, risk = run_case(
            try_name=args.try_name,
            r_value=r_value,
            terminals=terminals,
            params=params,
            reuse_existing=args.reuse_existing,
            debug=args.debug,
        )
        stage_log["cases"][str(r_value)] = {
            "backbone": backbone_summary,
            "risk_json": str(case_paths(args.try_name, r_value)["risk_json"]),
            "case_runtime_s": time.perf_counter() - case_t0,
        }
        rows.append(
            _case_row(
                r_request=r_value,
                label=f"OPT hierarchy R{r_value}",
                analysis=risk,
                original=original,
                backbone_summary=backbone_summary,
            )
        )

    stage_log["total_runtime_s"] = time.perf_counter() - t0
    manifest = write_summary(try_name=args.try_name, rows=rows, stage_log=stage_log, params=params)
    print(json.dumps({"rows": rows, "stage_log": stage_log, "manifest": str(manifest)}, indent=2))


if __name__ == "__main__":
    main()
