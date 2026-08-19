from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd

try:
    import contextily as ctx
except Exception:  # pragma: no cover - optional plotting dependency
    ctx = None

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.analyze_p2u_final_network import analyze_final_network  # noqa: E402
from figures.optimal_sfo.p2u_final_network import build_and_write_final_network  # noqa: E402
from figures.optimal_sfo.prepare_p2u_corridor_network import build_corridor_outputs  # noqa: E402
from figures.optimal_sfo.run_p2u_ilp_2edge import (  # noqa: E402
    SUPER_SOURCE,
    load_ilp_graph,
    solve_min_2edge,
    write_solution,
)
from figures.optimal_sfo.compare_p2u_old_new_reliability import (  # noqa: E402
    _analyze_one,
    _prepare_new_graph,
    _prepare_old_graph,
)


OUTPUT_DIR = ROOT / "outputs" / "optimal_hierarchy"
CORRIDOR_GPKG = ROOT / "outputs" / "optimal_sfo" / "p2u_terminal_corridors_road2_k10_3857.gpkg"


@dataclass(frozen=True)
class SweepParameters:
    p_mean: float = 5e-4
    ilp_time_limit_s: float = 600.0
    ilp_mip_gap: float = 0.05
    ilp_threads: int = 0
    ilp_max_cut_rounds: int = 100
    ilp_cut_mode: str = "callback"
    tree_mode: str = "street_forest"
    generalized_method: str = "projection"
    redundancy_mode: str = "exact"


def mode_output_dir(redundancy_mode: str) -> Path:
    return OUTPUT_DIR / ("exact" if redundancy_mode == "exact" else "max")


def manifest_path(redundancy_mode: str) -> Path:
    suffix = "exact" if redundancy_mode == "exact" else "max"
    return mode_output_dir(redundancy_mode) / f"p2u_hierarchical_redundancy_sweep_{suffix}.md"


def case_dir(r_value: int, redundancy_mode: str) -> Path:
    return mode_output_dir(redundancy_mode) / f"R{r_value:03d}"


def case_paths(r_value: int, redundancy_mode: str) -> dict[str, Path]:
    directory = case_dir(r_value, redundancy_mode)
    token = f"R{r_value}" if redundancy_mode == "exact" else f"Rmax{r_value}"
    return {
        "dir": directory,
        "backbone_gpkg": directory / f"p2u_backbone_{token}_3857.gpkg",
        "backbone_json": directory / f"p2u_backbone_{token}_summary.json",
        "backbone_md": directory / f"p2u_backbone_{token}_summary.md",
        "final_gpkg": directory / f"p2u_final_network_{token}_streetforest_3857.gpkg",
        "final_metadata": directory / f"p2u_final_network_{token}_streetforest_metadata.json",
        "topology_json": directory / f"p2u_final_network_{token}_topology.json",
        "topology_md": directory / f"p2u_final_network_{token}_topology.md",
        "risk_json": directory / f"p2u_final_network_{token}_risk.json",
        "network_png": directory / f"p2u_final_network_{token}_map.png",
    }


def _old_optimal_sfo_paths(r_max: int) -> dict[str, Path]:
    base = ROOT / "outputs" / "optimal_sfo"
    return {
        "backbone_gpkg": base / f"p2u_ilp_2edge_solution_Rmax{r_max}_3857.gpkg",
        "backbone_json": base / f"p2u_ilp_2edge_solution_Rmax{r_max}_summary.json",
        "backbone_md": base / f"p2u_ilp_2edge_solution_Rmax{r_max}_summary.md",
        "final_gpkg": base / f"p2u_final_network_Rmax{r_max}_streetforest_3857.gpkg",
        "final_metadata": base / f"p2u_final_network_Rmax{r_max}_streetforest_metadata.json",
    }


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _edge_key(u: str, v: str) -> tuple[str, str]:
    return tuple(sorted((str(u), str(v))))


def load_backbone_warm_start_edges(gpkg: Path | None) -> set[tuple[str, str]]:
    if gpkg is None or not gpkg.exists():
        return set()
    edges = gpd.read_file(gpkg, layer="solution_edges")
    warm_edges: set[tuple[str, str]] = set()
    for _, row in edges.iterrows():
        u = SUPER_SOURCE if str(row.terminal_kind_a) in {"source", "source_transformer"} else str(row.terminal_a)
        v = SUPER_SOURCE if str(row.terminal_kind_b) in {"source", "source_transformer"} else str(row.terminal_b)
        if u != v:
            warm_edges.add(_edge_key(u, v))
    return warm_edges


def prepare_corridors(*, rebuild: bool) -> dict[str, Any]:
    if CORRIDOR_GPKG.exists() and not rebuild:
        return {"status": "reused", "output_gpkg": str(CORRIDOR_GPKG)}
    t0 = time.perf_counter()
    summary = build_corridor_outputs()
    summary["stage_runtime_s"] = time.perf_counter() - t0
    return summary


def import_existing_optimal_sfo_case(r_max: int, redundancy_mode: str) -> dict[str, Any]:
    src = _old_optimal_sfo_paths(r_max)
    missing = [str(path) for path in src.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Cannot import R={r_max}; missing existing artifacts: {missing}")
    backbone = _read_json(src["backbone_json"])
    if redundancy_mode == "exact" and int(backbone.get("solution_cycle_rank", -1)) != int(r_max):
        raise ValueError(
            f"Cannot import {src['backbone_json']} as exact R={r_max}; "
            f"solution_cycle_rank={backbone.get('solution_cycle_rank')}"
        )

    dst = case_paths(r_max, redundancy_mode)
    dst["dir"].mkdir(parents=True, exist_ok=True)
    shutil.copy2(src["backbone_gpkg"], dst["backbone_gpkg"])
    shutil.copy2(src["backbone_md"], dst["backbone_md"])
    shutil.copy2(src["final_gpkg"], dst["final_gpkg"])

    backbone["output_gpkg"] = str(dst["backbone_gpkg"])
    backbone["imported_from"] = str(src["backbone_json"])
    backbone["redundancy_mode"] = redundancy_mode
    _write_json(dst["backbone_json"], backbone)

    metadata = _read_json(src["final_metadata"])
    metadata["input_backbone_gpkg"] = str(dst["backbone_gpkg"])
    metadata["output_gpkg"] = str(dst["final_gpkg"])
    metadata["metadata_json"] = str(dst["final_metadata"])
    metadata["imported_from"] = str(src["final_metadata"])
    _write_json(dst["final_metadata"], metadata)
    return {"status": "imported", "source": str(src["final_gpkg"])}


def solve_backbone(
    r_value: int,
    params: SweepParameters,
    *,
    reuse_existing: bool,
    warm_start_gpkg: Path | None = None,
) -> dict[str, Any]:
    paths = case_paths(r_value, params.redundancy_mode)
    if reuse_existing and paths["backbone_json"].exists() and paths["backbone_gpkg"].exists():
        return {**_read_json(paths["backbone_json"]), "stage_status": "reused"}

    graph, edges, transformer_nodes, source_nodes, original_source_incident = load_ilp_graph(CORRIDOR_GPKG)
    warm_start_edges = load_backbone_warm_start_edges(warm_start_gpkg)
    solution, summary = solve_min_2edge(
        graph,
        original_source_incident=original_source_incident,
        redundancy=r_value if params.redundancy_mode == "exact" else None,
        max_redundancy=r_value if params.redundancy_mode == "max" else None,
        time_limit=params.ilp_time_limit_s,
        mip_gap=params.ilp_mip_gap,
        threads=params.ilp_threads,
        max_cut_rounds=params.ilp_max_cut_rounds,
        cut_mode=params.ilp_cut_mode,
        warm_start_edges=warm_start_edges,
    )
    summary["warm_start_gpkg"] = str(warm_start_gpkg) if warm_start_gpkg else None
    summary["redundancy_mode"] = params.redundancy_mode
    write_solution(
        solution,
        summary,
        edges,
        transformer_nodes,
        source_nodes,
        output_gpkg=paths["backbone_gpkg"],
        summary_json=paths["backbone_json"],
        summary_md=paths["backbone_md"],
    )
    summary["stage_status"] = "solved"
    _write_json(paths["backbone_json"], summary)
    return summary


def build_final_network(r_value: int, params: SweepParameters, *, reuse_existing: bool) -> dict[str, Any]:
    paths = case_paths(r_value, params.redundancy_mode)
    if reuse_existing and paths["final_metadata"].exists() and paths["final_gpkg"].exists():
        return {**_read_json(paths["final_metadata"]), "stage_status": "reused"}
    t0 = time.perf_counter()
    metadata = build_and_write_final_network(
        backbone_summary=paths["backbone_json"],
        output_gpkg=paths["final_gpkg"],
        metadata_json=paths["final_metadata"],
        tree_mode=params.tree_mode,
    )
    metadata["stage_runtime_s"] = time.perf_counter() - t0
    metadata["stage_status"] = "built"
    _write_json(paths["final_metadata"], metadata)
    return metadata


def analyze_case(r_value: int, params: SweepParameters, *, reuse_existing: bool) -> dict[str, Any]:
    paths = case_paths(r_value, params.redundancy_mode)
    if reuse_existing and paths["risk_json"].exists() and paths["topology_json"].exists():
        return {**_read_json(paths["risk_json"]), "stage_status": "reused"}

    topology_t0 = time.perf_counter()
    topology = analyze_final_network(
        gpkg_path=paths["final_gpkg"],
        metadata_json=paths["final_metadata"],
        p_mean=params.p_mean,
        n_samples=0,
        output_json=paths["topology_json"],
        output_md=paths["topology_md"],
    )
    topology["stage_runtime_s"] = time.perf_counter() - topology_t0
    _write_json(paths["topology_json"], topology)

    graph, sources, extra = _prepare_new_graph(paths["final_gpkg"], paths["final_metadata"])
    risk = _analyze_one(
        name=f"p2u_hierarchical_{params.redundancy_mode}_R{r_value}",
        graph=graph,
        sources=sources,
        extra=extra,
        p_mean=params.p_mean,
        generalized=True,
        generalized_method=params.generalized_method,
    )
    risk["topology_summary"] = topology
    risk["stage_status"] = "computed"
    _write_json(paths["risk_json"], risk)
    return risk


def compute_original_reference(params: SweepParameters, *, reuse_existing: bool) -> dict[str, Any]:
    output = OUTPUT_DIR / "p2u_original_reference_risk.json"
    if reuse_existing and output.exists():
        return _read_json(output)
    graph, sources, extra = _prepare_old_graph()
    reference = _analyze_one(
        name="original_full_p2u_mv",
        graph=graph,
        sources=sources,
        extra=extra,
        p_mean=params.p_mean,
        generalized=True,
        generalized_method=params.generalized_method,
    )
    _write_json(output, reference)
    return reference


def _risk_value(analysis: dict[str, Any], key: str) -> float:
    return float(analysis["risk_float"].get(key, 0.0))


def _topology_float(topology: dict[str, Any], key: str) -> float | None:
    value = topology.get(key)
    if value is None:
        return None
    return float(value)


def _case_row(
    *,
    r_request: int | None,
    label: str,
    analysis: dict[str, Any],
    original: dict[str, Any],
    backbone_summary: dict[str, Any] | None,
    is_original: bool = False,
) -> dict[str, Any]:
    topology = analysis["topology"]
    topology_summary = analysis.get("topology_summary", {})
    risk_total = _risk_value(analysis, "total")
    op = _risk_value(analysis, "tree") + _risk_value(analysis, "nonbridge_section")
    op2 = _risk_value(analysis, "internal_regular_chains") + _risk_value(analysis, "structural_generalized_chains")

    if is_original:
        length_m = float(topology["total_length_m"])
    else:
        length_m = float(topology_summary.get("total_physical_length_m", topology.get("total_length_m", 0.0)))

    original_length = float(original["topology"]["total_length_m"])
    original_r = float(original["topology"]["source_contracted_cycle_rank_R"])
    original_risk = _risk_value(original, "total")
    original_op = _risk_value(original, "tree") + _risk_value(original, "nonbridge_section")
    original_op2 = _risk_value(original, "internal_regular_chains") + _risk_value(original, "structural_generalized_chains")

    achieved_r = float(topology["source_contracted_cycle_rank_R"])
    return {
        "network": label,
        "r_request": r_request,
        "length_km": length_m / 1000.0,
        "r_theory": achieved_r,
        "z_w": 1.0 - length_m / original_length,
        "z_r": 1.0 - achieved_r / original_r,
        "z_f": 1.0 - risk_total / original_risk,
        "z_f_p": 1.0 - op / original_op if original_op else None,
        "z_f_p2": 1.0 - op2 / original_op2 if original_op2 else None,
        "risk_total": risk_total,
        "risk_o_p": op,
        "risk_o_p2": op2,
        "risk_tree": _risk_value(analysis, "tree"),
        "risk_section": _risk_value(analysis, "nonbridge_section"),
        "risk_internal": _risk_value(analysis, "internal_regular_chains"),
        "risk_structural": _risk_value(analysis, "structural_generalized_chains"),
        "bridges": int(topology["source_contracted_bridges"]),
        "structure_nodes": int(topology["structure_graph_nodes"]),
        "regular_chains": int(topology["regular_chains"]),
        "generalized_chains": int(topology["generalized_chains"]),
        "gen_lambda_mean_km": (
            None
            if _topology_float(topology, "generalized_chain_mean_effective_lambda_m") is None
            else _topology_float(topology, "generalized_chain_mean_effective_lambda_m") / 1000.0
        ),
        "gen_lambda_sigma_over_mean": _topology_float(
            topology, "generalized_chain_std_over_mean_effective_lambda"
        ),
        "gen_lambda_max_km": (
            None
            if _topology_float(topology, "generalized_chain_max_effective_lambda_m") is None
            else _topology_float(topology, "generalized_chain_max_effective_lambda_m") / 1000.0
        ),
        "gen_physical_lambda_mean_km": (
            None
            if _topology_float(topology, "generalized_chain_mean_length_m") is None
            else _topology_float(topology, "generalized_chain_mean_length_m") / 1000.0
        ),
        "gen_physical_lambda_sigma_over_mean": _topology_float(
            topology, "generalized_chain_std_over_mean_length"
        ),
        "gen_physical_lambda_max_km": (
            None
            if _topology_float(topology, "generalized_chain_max_length_m") is None
            else _topology_float(topology, "generalized_chain_max_length_m") / 1000.0
        ),
        "ilp_status": None if backbone_summary is None else backbone_summary.get("status_name"),
        "ilp_runtime_s": None if backbone_summary is None else backbone_summary.get("runtime_s"),
        "ilp_mip_gap": None if backbone_summary is None else backbone_summary.get("mip_gap"),
        "decomposition_runtime_s": analysis.get("decomposition_runtime_seconds"),
    }


def classify_backbone_edges(backbone: gpd.GeoDataFrame) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    graph = nx.Graph()
    for idx, row in backbone.reset_index(drop=True).iterrows():
        u = str(row.terminal_a)
        v = str(row.terminal_b)
        if u != v:
            graph.add_edge(u, v, weight=float(row.length_m), row_id=int(idx))
    closed_ids: set[int] = set()
    for component in nx.connected_components(graph):
        sub = graph.subgraph(component)
        for _, _, data in nx.minimum_spanning_edges(sub, data=True):
            closed_ids.add(int(data["row_id"]))
    closed = backbone.reset_index(drop=True).loc[lambda df: df.index.isin(closed_ids)].copy()
    open_ties = backbone.reset_index(drop=True).loc[lambda df: ~df.index.isin(closed_ids)].copy()
    return closed, open_ties


def plot_network_png(r_value: int, params: SweepParameters, *, reuse_existing: bool) -> dict[str, Any]:
    paths = case_paths(r_value, params.redundancy_mode)
    if reuse_existing and paths["network_png"].exists():
        return {"stage_status": "reused", "output_png": str(paths["network_png"])}

    t0 = time.perf_counter()
    backbone = gpd.read_file(paths["final_gpkg"], layer="final_backbone_edges").to_crs("EPSG:3857")
    tree = gpd.read_file(paths["final_gpkg"], layer="final_tree_attachment_edges").to_crs("EPSG:3857")
    transformers = gpd.read_file(paths["final_gpkg"], layer="final_transformer_nodes").to_crs("EPSG:3857")
    sources = gpd.read_file(paths["final_gpkg"], layer="final_source_nodes").to_crs("EPSG:3857")
    closed, open_ties = classify_backbone_edges(backbone)

    fig, ax = plt.subplots(figsize=(7.0, 7.0), dpi=240)
    if len(tree):
        tree.plot(ax=ax, color="#888888", linewidth=0.28, alpha=0.55, zorder=3)
    if len(closed):
        closed.plot(ax=ax, color="#bf3b30", linewidth=0.55, alpha=0.95, zorder=4)
    if len(open_ties):
        open_ties.plot(ax=ax, color="#f2b21b", linewidth=1.0, linestyle=(0, (3, 2)), alpha=1.0, zorder=5)
    if len(transformers):
        transformers.plot(ax=ax, facecolor="none", edgecolor="#222222", markersize=1.2, linewidth=0.25, alpha=0.75, zorder=6)
    if len(sources):
        sources.plot(ax=ax, color="#bf3b30", edgecolor="white", markersize=24, linewidth=0.45, zorder=7)

    bounds = pd.concat([backbone.geometry, tree.geometry], ignore_index=True).total_bounds
    pad_x = (bounds[2] - bounds[0]) * 0.04
    pad_y = (bounds[3] - bounds[1]) * 0.04
    ax.set_xlim(bounds[0] - pad_x, bounds[2] + pad_x)
    ax.set_ylim(bounds[1] - pad_y, bounds[3] + pad_y)
    if ctx is not None:
        try:
            ctx.add_basemap(ax, crs="EPSG:3857", source=ctx.providers.CartoDB.Positron, attribution=False)
        except Exception as exc:  # pragma: no cover - depends on tile/network availability
            ax.text(0.01, 0.01, f"Basemap unavailable: {exc.__class__.__name__}", transform=ax.transAxes, fontsize=6)
    ax.set_axis_off()
    symbol = "=" if params.redundancy_mode == "exact" else "<="
    ax.set_title(f"P2U hierarchical network, R {symbol} {r_value}", fontsize=10)
    paths["network_png"].parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(paths["network_png"], bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return {
        "stage_status": "plotted",
        "output_png": str(paths["network_png"]),
        "plot_runtime_s": time.perf_counter() - t0,
        "normally_closed_backbone_edges": int(len(closed)),
        "normally_open_tie_edges": int(len(open_ties)),
    }


def write_summary(rows: list[dict[str, Any]], stage_log: dict[str, Any], params: SweepParameters) -> None:
    out_dir = mode_output_dir(params.redundancy_mode)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    table_csv = out_dir / f"p2u_hierarchical_redundancy_sweep_{params.redundancy_mode}_table.csv"
    table_json = out_dir / f"p2u_hierarchical_redundancy_sweep_{params.redundancy_mode}_table.json"
    table_csv.write_text(df.to_csv(index=False, lineterminator="\n"), encoding="utf-8", newline="\n")
    _write_json(table_json, {"rows": rows, "stage_log": stage_log})
    diagnostic_png = write_diagnostic_plot(df, params)

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
    md_table = _markdown_table(df[display_cols])
    r_values = [str(int(row["r_request"])) for row in rows if row.get("r_request") is not None]
    achieved_r_values = sorted({int(row["r_theory"]) for row in rows if row.get("r_request") is not None})
    flat_budget_note = ""
    if params.redundancy_mode == "max" and len(achieved_r_values) == 1 and len(r_values) > 1:
        flat_budget_note = (
            f"- All requested upper bounds selected the same achieved `R={achieved_r_values[0]}` backbone. "
            "This is expected for a pure minimum-length objective with `R <= R_max`: extra redundancy is allowed but not rewarded."
        )
    exact_note = ""
    if params.redundancy_mode == "exact":
        exact_note = "- This run enforces exact source-contracted redundancy, so `r_request` should match `r_theory` for successful rows."
    reproduce_command = (
        "& C:\\Users\\rotem\\anaconda3\\envs\\reliability\\python.exe "
        f"figures\\optimal_hierarchy\\redundancy_sweep.py --redundancy-mode {params.redundancy_mode} --r-values {' '.join(r_values or ['50'])} --reuse-existing"
    )
    constraint_text = (
        "$$" + "\n" + "|E_\\mathrm{selected}| - |V| + 1 = R_\\mathrm{target}." + "\n" + "$$"
        if params.redundancy_mode == "exact"
        else "$$" + "\n" + "|E_\\mathrm{selected}| - |V| + 1 \\le R_\\max." + "\n" + "$$"
    )
    constraint_label = "target" if params.redundancy_mode == "exact" else "upper-bound"
    lines = [
        f"# P2U Hierarchical Redundancy Sweep ({params.redundancy_mode})",
        "",
        "## Results",
        "",
        f"**Goal**: prepare an article-level P2U experiment where a road-constrained 2-edge-connected backbone is built under a redundancy {constraint_label} and all remaining transformers are attached as a street forest.",
        "",
        "**Main result**: this manifest is the first clean sweep scaffold. It uses the R50-style parameters that produced the high `O(p^2)` structural-risk example and records the achieved network for each requested redundancy budget.",
        "",
        "All reliability rows use deterministic decomposition with length-scaled `p_mean = 5e-4`.",
        "",
        md_table,
        "",
        "**Figures**:",
        "",
    ]
    for row in rows:
        if row["network"] == "Original P2U MV":
            continue
        png = case_paths(int(row["r_request"]), params.redundancy_mode)["network_png"]
        if png.exists():
            rel = png.relative_to(out_dir)
            lines.append(f"- `{row['network']}`: ![]({rel.as_posix()})")
    lines.extend(
        [
            "",
            "**Insights**:",
            "",
            "- `O(p) = Tree + Section` is the first-order bridge/section risk.",
            "- `O(p^2) = Internal + Structural` is the second-order chain and generalized-chain risk.",
            "- `gen_lambda_*` reports demand-normalized generalized-chain effective lengths, using `tilde_lambda_q = lambda_q sqrt(Q w_q / W)`.",
            "- The R50 case is kept as the first benchmark because it shows that a low-redundancy 2-connected backbone can leave a large `O(p^2)` component.",
            exact_note,
            flat_budget_note,
            "",
            "**What this does not show**:",
            "",
            "- This stage produces a simulation table and preview maps, not the final article figure layout.",
            "- ILP results are only as strong as the recorded Gurobi status and MIP gap for each row.",
            "- The experiment remains connectivity reliability, not voltage/power-flow validation.",
            "",
            "**Reproduce**:",
            "",
            "```powershell",
            reproduce_command,
            "```",
            "",
            "## Algorithm",
            "",
            "For each requested redundancy value, the pipeline uses the road-corridor terminal graph and solves a minimum-length 2-edge-connected backbone ILP:",
            "",
            "$$",
            "\\min \\sum_e w_e x_e",
            "$$",
            "",
            "subject to degree, source-incidence, lazy 2-edge cut constraints, and the selected redundancy constraint:",
            "",
            constraint_text,
            "",
            "The final graph attaches non-backbone transformer terminals using shortest paths on the street network and contracts street/transformer chains into switch-section edges.",
            "",
            "The reliability split is:",
            "",
            "$$",
            "F = F_{O(p)} + F_{O(p^2)}, \\quad F_{O(p)} = F_\\mathrm{tree} + F_\\mathrm{section}, \\quad F_{O(p^2)} = F_\\mathrm{internal} + F_\\mathrm{structural}.",
            "$$",
            "",
            "Relative indexes are computed against the original P2U MV network:",
            "",
            "$$",
            "Z_W = 1 - W/W_0, \\quad Z_R = 1 - R/R_0, \\quad Z_F = 1 - F/F_0.",
            "$$",
            "",
            "Current parameters:",
            "",
            f"- `p_mean = {params.p_mean}`",
            f"- `time_limit = {params.ilp_time_limit_s}` s",
            f"- `MIPGap = {params.ilp_mip_gap}`",
            f"- `Threads = {params.ilp_threads}`",
            f"- `cut_mode = {params.ilp_cut_mode}`",
            f"- `tree_mode = {params.tree_mode}`",
            f"- `generalized_method = {params.generalized_method}`",
            f"- `redundancy_mode = {params.redundancy_mode}`",
            "",
            "## Implementation",
            "",
            "Entry point:",
            "",
            "- `figures/optimal_hierarchy/redundancy_sweep.py`",
            "",
            "Shared code reused from the exploratory P2U implementation:",
            "",
            "- `figures/optimal_sfo/prepare_p2u_corridor_network.py` builds the road-corridor candidate graph.",
            "- `figures/optimal_sfo/run_p2u_ilp_2edge.py` solves the backbone ILP.",
            "- `figures/optimal_sfo/p2u_final_network.py` builds the final backbone-plus-forest graph.",
            "- `figures/optimal_sfo/analyze_p2u_final_network.py` computes topology, load, and length summaries.",
            "- `figures/optimal_sfo/compare_p2u_old_new_reliability.py` provides the deterministic risk-decomposition path.",
            "",
            "Outputs:",
            "",
            f"- `{table_csv.relative_to(ROOT)}`",
            f"- `{table_json.relative_to(ROOT)}`",
            f"- `{manifest_path(params.redundancy_mode).relative_to(ROOT)}`",
            f"- `{diagnostic_png.relative_to(ROOT)}`",
            f"- `{out_dir.relative_to(ROOT)}/R*/p2u_final_network_*_map.png`",
            "",
            "Stage log:",
            "",
            "```json",
            json.dumps(stage_log, indent=2),
            "```",
            "",
            "No project research log was found at `codex/documents/research_logs/research_log.md`, so no research-log entry was updated.",
            "",
        ]
    )
    manifest_path(params.redundancy_mode).write_text("\n".join(lines), encoding="utf-8")


def write_diagnostic_plot(df: pd.DataFrame, params: SweepParameters) -> Path:
    out_dir = mode_output_dir(params.redundancy_mode)
    out_png = out_dir / f"p2u_hierarchical_redundancy_sweep_{params.redundancy_mode}_diagnostic.png"
    sweep = df[df["r_request"].notna()].copy()
    if sweep.empty:
        return out_png
    sweep["r_request"] = sweep["r_request"].astype(float)
    sweep = sweep.sort_values("r_request")

    panels = [
        ("z_w", "$Z_W$"),
        ("z_f", "$Z_F$"),
        ("z_f_p", "$Z_{F,p}$"),
        ("z_f_p2", "$Z_{F,p^2}$"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.0), dpi=220, sharex=True)
    for ax, (column, title) in zip(axes.ravel(), panels):
        ax.plot(
            sweep["r_request"],
            sweep[column],
            color="#222222",
            marker="o",
            markersize=3.5,
            linewidth=1.2,
        )
        ax.axhline(0.0, color="#999999", linewidth=0.7, linestyle=":")
        ax.set_title(title, fontsize=10)
        ax.grid(True, color="#dddddd", linewidth=0.5, alpha=0.8)
        ax.tick_params(labelsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel("R", fontsize=9)
    fig.suptitle(f"P2U exact-R hierarchy sweep ({params.redundancy_mode})", fontsize=11)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    return out_png


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
    parser = argparse.ArgumentParser(description="Run the article-level P2U hierarchical redundancy sweep.")
    parser.add_argument("--r-values", type=int, nargs="+", default=[50])
    parser.add_argument("--rebuild-corridors", action="store_true")
    parser.add_argument("--reuse-existing", action="store_true")
    parser.add_argument(
        "--import-existing-optimal-sfo",
        action="store_true",
        help="Import matching existing exploratory optimal_sfo artifacts before analysis; useful for bootstrapping R50/R171 without rerunning ILP.",
    )
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--mip-gap", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--max-cut-rounds", type=int, default=100)
    parser.add_argument("--cut-mode", choices=["iterative", "callback"], default="callback")
    parser.add_argument("--p-mean", type=float, default=5e-4)
    parser.add_argument("--generalized-method", choices=["projection", "networkx"], default="projection")
    parser.add_argument("--redundancy-mode", choices=["exact", "max"], default="exact")
    parser.add_argument(
        "--refresh-analysis",
        action="store_true",
        help="Reuse existing corridor/backbone/final-network artifacts, but recompute reliability/topology analysis JSON and tables.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = SweepParameters(
        p_mean=args.p_mean,
        ilp_time_limit_s=args.time_limit,
        ilp_mip_gap=args.mip_gap,
        ilp_threads=args.threads,
        ilp_max_cut_rounds=args.max_cut_rounds,
        ilp_cut_mode=args.cut_mode,
        generalized_method=args.generalized_method,
        redundancy_mode=args.redundancy_mode,
    )
    mode_output_dir(params.redundancy_mode).mkdir(parents=True, exist_ok=True)
    stage_log: dict[str, Any] = {"parameters": params.__dict__, "cases": {}}
    t0 = time.perf_counter()

    stage_log["corridors"] = prepare_corridors(rebuild=args.rebuild_corridors)
    analysis_reuse_existing = args.reuse_existing and not args.refresh_analysis
    original = compute_original_reference(params, reuse_existing=analysis_reuse_existing)
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

    previous_backbone_gpkg: Path | None = None
    for r_max in args.r_values:
        paths = case_paths(r_max, params.redundancy_mode)
        paths["dir"].mkdir(parents=True, exist_ok=True)
        case_t0 = time.perf_counter()
        log: dict[str, Any] = {}
        if args.import_existing_optimal_sfo:
            log["import"] = import_existing_optimal_sfo_case(r_max, params.redundancy_mode)
        backbone = solve_backbone(
            r_max,
            params,
            reuse_existing=args.reuse_existing or args.import_existing_optimal_sfo,
            warm_start_gpkg=previous_backbone_gpkg,
        )
        if paths["backbone_gpkg"].exists():
            previous_backbone_gpkg = paths["backbone_gpkg"]
        else:
            log["backbone"] = backbone
            log["case_runtime_s"] = time.perf_counter() - case_t0
            stage_log["cases"][str(r_max)] = log
            raise RuntimeError(f"No backbone GeoPackage was written for R={r_max}: {backbone}")
        final = build_final_network(r_max, params, reuse_existing=args.reuse_existing or args.import_existing_optimal_sfo)
        risk = analyze_case(r_max, params, reuse_existing=analysis_reuse_existing)
        plot = None if args.skip_plots else plot_network_png(r_max, params, reuse_existing=args.reuse_existing)
        log.update(
            {
                "backbone": backbone,
                "final_network": final,
                "risk_json": str(paths["risk_json"]),
                "plot": plot,
                "case_runtime_s": time.perf_counter() - case_t0,
            }
        )
        stage_log["cases"][str(r_max)] = log
        rows.append(
            _case_row(
                r_request=r_max,
                label=f"Hierarchical road R{r_max}",
                analysis=risk,
                original=original,
                backbone_summary=backbone,
            )
        )

    stage_log["total_runtime_s"] = time.perf_counter() - t0
    write_summary(rows, stage_log, params)
    print(json.dumps({"rows": rows, "stage_log": stage_log, "manifest": str(manifest_path(params.redundancy_mode))}, indent=2))


if __name__ == "__main__":
    main()
