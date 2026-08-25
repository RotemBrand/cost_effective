from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from shapely.geometry import Point

try:
    import contextily as ctx
except Exception:  # pragma: no cover - optional map tiles
    ctx = None

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_hierarchy.redundancy_sweep import case_paths  # noqa: E402
from figures.optimal_sfo.compare_p2u_old_new_reliability import (  # noqa: E402
    _prepare_new_graph,
    _prepare_old_graph,
)
from figures.p2u_full_mv_reliability import build_full_mv_graph  # noqa: E402
from indexes import GraphRel  # noqa: E402


OUTPUT_DIR = ROOT / "outputs" / "optimal_hierarchy" / "exact"

TREE_COLOR = "#2f855a"
BACKBONE_COLOR = "#bf3b30"
TIE_COLOR = "#f2b21b"
CHAIN_BASE_COLOR = "#9ca3af"


@dataclass(frozen=True)
class NetworkCase:
    key: str
    label: str
    graph: nx.Graph
    sources: list[Any]
    gpkg: Path | None = None


def _terminal_edge_key(u: Any, v: Any) -> tuple[str, str]:
    return tuple(sorted((str(u).strip(), str(v).strip()), key=repr))


def _source_contracted_edge_classes(graph: nx.Graph, sources: list[Any]) -> tuple[set[tuple[Any, Any]], set[tuple[Any, Any]]]:
    source_set = set(sources)
    contracted = nx.Graph()
    edge_map: dict[tuple[Any, Any], list[tuple[Any, Any]]] = {}
    for u, v in graph.edges:
        cu = "__SOURCE__" if u in source_set else u
        cv = "__SOURCE__" if v in source_set else v
        if cu == cv:
            continue
        key = tuple(sorted((cu, cv), key=repr))
        contracted.add_edge(*key)
        edge_map.setdefault(key, []).append(tuple(sorted((u, v), key=repr)))
    bridge_keys = {tuple(sorted(edge, key=repr)) for edge in nx.bridges(contracted)}
    tree_edges = {edge for key in bridge_keys for edge in edge_map.get(key, [])}
    backbone_edges = {edge for key, originals in edge_map.items() if key not in bridge_keys for edge in originals}
    return tree_edges, backbone_edges


def _original_snapshot_tables() -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, NetworkCase]:
    _raw_graph, mv_rows, _switches = build_full_mv_graph()
    old_graph, sources, _extra = _prepare_old_graph()
    tree_edges, backbone_edges = _source_contracted_edge_classes(old_graph, sources)
    rows = mv_rows.copy()
    rows["_edge_key"] = [
        tuple(sorted((str(row.NodeA).strip(), str(row.NodeB).strip()), key=repr))
        for _, row in rows.iterrows()
    ]
    rows["is_tie"] = rows["Status"].astype(str).str.strip() == "0"
    rows["component_class"] = np.where(rows["_edge_key"].isin(tree_edges), "tree", "backbone")
    tree = rows[(rows["component_class"] == "tree") & (~rows["is_tie"])].to_crs("EPSG:3857")
    backbone = rows[(rows["component_class"] == "backbone") & (~rows["is_tie"])].to_crs("EPSG:3857")
    ties = rows[rows["is_tie"]].to_crs("EPSG:3857")

    node_rows = []
    seen = set()
    for _, row in rows.iterrows():
        coords = list(row.geometry.coords)
        for node, xy in ((str(row.NodeA).strip(), coords[0]), (str(row.NodeB).strip(), coords[-1])):
            if node in seen:
                continue
            seen.add(node)
            data = old_graph.nodes.get(node, {})
            node_rows.append(
                {
                    "terminal_id": node,
                    "size_kva": float(data.get("weight", 0.0)),
                    "is_source": bool(data.get("is_source", False)),
                    "geometry": Point(xy),
                }
            )
    nodes = gpd.GeoDataFrame(node_rows, geometry="geometry", crs=mv_rows.crs).to_crs("EPSG:3857")
    source_ids = set(sources)
    nodes["is_source"] = nodes["terminal_id"].isin(source_ids)
    return tree, backbone, ties, nodes, NetworkCase("original", "Original", old_graph, sources, None)


def _minimum_tree_edge_indices(edges: gpd.GeoDataFrame) -> set[int]:
    graph = nx.Graph()
    for idx, row in edges.reset_index(drop=True).iterrows():
        u = str(row.terminal_a)
        v = str(row.terminal_b)
        if u != v:
            graph.add_edge(u, v, weight=float(row.length_m), row_id=int(idx))
    selected: set[int] = set()
    for component in nx.connected_components(graph):
        sub = graph.subgraph(component)
        for _, _, data in nx.minimum_spanning_edges(sub, data=True):
            selected.add(int(data["row_id"]))
    return selected


def _final_snapshot_tables(
    gpkg: Path,
    metadata: Path,
    key: str,
    label: str,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, NetworkCase]:
    backbone = gpd.read_file(gpkg, layer="final_backbone_edges").to_crs("EPSG:3857")
    tree = gpd.read_file(gpkg, layer="final_tree_attachment_edges").to_crs("EPSG:3857")
    transformers = gpd.read_file(gpkg, layer="final_transformer_nodes").to_crs("EPSG:3857")
    sources = gpd.read_file(gpkg, layer="final_source_nodes").to_crs("EPSG:3857")
    graph, source_ids, _ = _prepare_new_graph(gpkg, metadata)
    closed_ids = _minimum_tree_edge_indices(backbone)
    ties = backbone.loc[lambda df: ~df.index.isin(closed_ids)].copy()
    closed = backbone.loc[lambda df: df.index.isin(closed_ids)].copy()
    nodes = pd.concat([transformers, sources], ignore_index=True)
    nodes["is_source"] = nodes.get("kind", "").astype(str).eq("source")
    return tree, closed, ties, gpd.GeoDataFrame(nodes, geometry="geometry", crs=transformers.crs), NetworkCase(key, label, graph, source_ids, gpkg)


def _decompose(case: NetworkCase):
    node_weights = {node: float(data.get("weight", 0.0)) for node, data in case.graph.nodes(data=True)}
    return GraphRel(case.graph, nodes_weight=node_weights, sources=case.sources).decompose(
        include_generalized_chains=True,
        generalized_component_method="projection",
    )


def _snapshot_edge_geometry_table(
    tree: gpd.GeoDataFrame,
    backbone: gpd.GeoDataFrame,
    ties: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    tables = []
    for source in (tree, backbone, ties):
        if source.empty:
            continue
        table = source.copy()
        if {"terminal_a", "terminal_b"}.issubset(table.columns):
            table["endpoint_a"] = table["terminal_a"].astype(str)
            table["endpoint_b"] = table["terminal_b"].astype(str)
        elif {"NodeA", "NodeB"}.issubset(table.columns):
            table["endpoint_a"] = table["NodeA"].astype(str).str.strip()
            table["endpoint_b"] = table["NodeB"].astype(str).str.strip()
        else:
            continue
        tables.append(table[["endpoint_a", "endpoint_b", "geometry"]])
    if not tables:
        return gpd.GeoDataFrame(columns=["endpoint_a", "endpoint_b", "geometry"], geometry="geometry", crs="EPSG:3857")
    return gpd.GeoDataFrame(pd.concat(tables, ignore_index=True), geometry="geometry", crs=tables[0].crs)


def _chain_effective_lambda_km(chain: Any, *, q: int, total_weight: float) -> float:
    if q <= 0 or total_weight <= 0.0:
        return 0.0
    return float(chain.length) * math.sqrt(q * max(float(chain.total_weight), 0.0) / total_weight) / 1000.0


def generalized_chain_edge_table(
    network_case: NetworkCase,
    tree: gpd.GeoDataFrame,
    backbone: gpd.GeoDataFrame,
    ties: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    edge_rows = _snapshot_edge_geometry_table(tree, backbone, ties)
    if edge_rows.empty:
        return edge_rows.assign(chain_id=pd.Series(dtype=int), lambda_km=pd.Series(dtype=float))
    geometry_by_key = {
        _terminal_edge_key(row.endpoint_a, row.endpoint_b): row.geometry
        for _, row in edge_rows.iterrows()
    }
    decomp = _decompose(network_case)
    source_edge_to_original = {
        tuple(sorted((u, v), key=repr)): _terminal_edge_key(*data["original_edge"])
        for u, v, data in decomp.source_graph.edges(data=True)
        if "original_edge" in data
    }
    rows = []
    q = len(decomp.generalized_chains)
    for chain_id, chain in enumerate(decomp.generalized_chains):
        seen: set[tuple[str, str]] = set()
        lambda_km = _chain_effective_lambda_km(chain, q=q, total_weight=float(decomp.total_weight))
        for source_edge in chain.original_edges:
            terminal_key = source_edge_to_original.get(tuple(sorted(source_edge, key=repr)))
            if terminal_key is None or terminal_key in seen:
                continue
            geometry = geometry_by_key.get(terminal_key)
            if geometry is None:
                continue
            seen.add(terminal_key)
            rows.append(
                {
                    "chain_id": int(chain_id),
                    "lambda_km": float(lambda_km),
                    "endpoint_a": terminal_key[0],
                    "endpoint_b": terminal_key[1],
                    "geometry": geometry,
                }
            )
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=edge_rows.crs)


def chain_width_scale(chain_tables: list[gpd.GeoDataFrame], *, quantile: float = 0.98) -> dict[str, float]:
    values = []
    for table in chain_tables:
        if table.empty or "lambda_km" not in table.columns:
            continue
        chain_values = table[["chain_id", "lambda_km"]].drop_duplicates()["lambda_km"].astype(float)
        values.extend(float(value) for value in chain_values if math.isfinite(float(value)))
    if not values:
        return {"clip_km": 1.0, "min_width": 0.36, "max_width": 2.35}
    clip = float(np.quantile(np.asarray(values, dtype=float), quantile))
    return {"clip_km": max(clip, 1e-9), "min_width": 0.36, "max_width": 2.35}


def _chain_line_width(lambda_km: float, scale: dict[str, float]) -> float:
    clipped = min(max(float(lambda_km), 0.0), float(scale["clip_km"]))
    fraction = clipped / float(scale["clip_km"]) if float(scale["clip_km"]) > 0.0 else 0.0
    return float(scale["min_width"]) + fraction * (float(scale["max_width"]) - float(scale["min_width"]))


def _plot_sources(ax: plt.Axes, sources: gpd.GeoDataFrame) -> None:
    if len(sources):
        sources.plot(ax=ax, color=BACKBONE_COLOR, edgecolor="white", markersize=9, linewidth=0.35, zorder=7)


def _plot_component_snapshot(
    ax: plt.Axes,
    tree: gpd.GeoDataFrame,
    backbone: gpd.GeoDataFrame,
    ties: gpd.GeoDataFrame,
    nodes: gpd.GeoDataFrame,
) -> None:
    if len(tree):
        tree.plot(ax=ax, color=TREE_COLOR, lw=0.32, alpha=0.62, zorder=2)
    if len(backbone):
        backbone.plot(ax=ax, color=BACKBONE_COLOR, lw=0.58, alpha=0.98, zorder=3)
    _plot_sources(ax, nodes[nodes["is_source"].astype(bool)].copy())
    if len(ties):
        ties.plot(ax=ax, color=TIE_COLOR, lw=1.35, linestyle=(0, (3, 2)), alpha=1.0, zorder=8)


def _plot_generalized_chain_snapshot(
    ax: plt.Axes,
    backbone: gpd.GeoDataFrame,
    ties: gpd.GeoDataFrame,
    chain_table: gpd.GeoDataFrame,
    *,
    width_scale: dict[str, float],
) -> None:
    if len(backbone):
        backbone.plot(ax=ax, color=CHAIN_BASE_COLOR, lw=0.35, alpha=0.42, zorder=2)
    if len(ties):
        ties.plot(ax=ax, color=CHAIN_BASE_COLOR, lw=0.35, alpha=0.42, zorder=2)
    if chain_table.empty:
        return
    cmap = plt.get_cmap("tab20")
    for chain_id, group in chain_table.groupby("chain_id", sort=True):
        color = mpl.colors.to_hex(cmap(int(chain_id) % cmap.N))
        lambda_km = float(group["lambda_km"].iloc[0])
        group.plot(ax=ax, color=color, lw=_chain_line_width(lambda_km, width_scale), alpha=0.98, zorder=4)


def plot_cases(
    cases: list[tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, NetworkCase]],
    *,
    output_prefix: Path,
) -> dict[str, Any]:
    fig, axes = plt.subplots(2, len(cases), figsize=(3.6 * len(cases), 4.45), dpi=220)
    if len(cases) == 1:
        axes = np.asarray(axes).reshape(2, 1)
    chain_tables = [
        generalized_chain_edge_table(network_case, tree, backbone, ties)
        for tree, backbone, ties, _nodes, network_case in cases
    ]
    width_scale = chain_width_scale(chain_tables)
    bounds_list = []
    sidecar_cases = []
    for col, (tree, backbone, ties, nodes, network_case) in enumerate(cases):
        _plot_component_snapshot(axes[0, col], tree, backbone, ties, nodes)
        _plot_generalized_chain_snapshot(axes[1, col], backbone, ties, chain_tables[col], width_scale=width_scale)
        axes[0, col].set_title(network_case.label, fontsize=10)
        axes[1, col].set_title("2-component generalized chains", fontsize=8.5)
        axes[0, col].set_axis_off()
        axes[1, col].set_axis_off()
        geoms = [g.geometry for g in (tree, backbone, ties) if len(g)]
        if geoms:
            bounds_list.append(pd.concat(geoms, ignore_index=True).total_bounds)
        lambda_values = chain_tables[col][["chain_id", "lambda_km"]].drop_duplicates()["lambda_km"].astype(float) if not chain_tables[col].empty else pd.Series(dtype=float)
        sidecar_cases.append(
            {
                "key": network_case.key,
                "label": network_case.label,
                "generalized_chain_count": int(len(lambda_values)),
                "lambda_mean_km": None if lambda_values.empty else float(lambda_values.mean()),
                "lambda_max_km": None if lambda_values.empty else float(lambda_values.max()),
            }
        )

    if bounds_list:
        bounds = np.asarray(bounds_list)
        xmin, ymin = np.nanmin(bounds[:, :2], axis=0)
        xmax, ymax = np.nanmax(bounds[:, 2:], axis=0)
        pad_x = (xmax - xmin) * 0.03
        pad_y = (ymax - ymin) * 0.03
        for ax in axes.ravel():
            ax.set_xlim(xmin - pad_x, xmax + pad_x)
            ax.set_ylim(ymin - pad_y, ymax + pad_y)
            if ctx is not None:
                try:
                    ctx.add_basemap(ax, crs="EPSG:3857", source=ctx.providers.CartoDB.Positron, attribution=False)
                except Exception:
                    pass

    handles = [
        Line2D([0], [0], color=BACKBONE_COLOR, lw=1.5, label="2-component/backbone"),
        Line2D([0], [0], color=TREE_COLOR, lw=1.5, label="1-component/tree"),
        Line2D([0], [0], color=TIE_COLOR, lw=2.0, ls="--", label="normally open ties"),
        Line2D([0], [0], color=CHAIN_BASE_COLOR, lw=1.3, label="2-component base"),
        Line2D([0], [0], color=plt.get_cmap("tab20")(0), lw=1.8, label=r"generalized chain, width $\propto \tilde{\lambda}$"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.01))
    fig.subplots_adjust(bottom=0.16, top=0.92, wspace=0.02, hspace=0.10)
    png = output_prefix.with_suffix(".png")
    svg = output_prefix.with_suffix(".svg")
    sidecar = output_prefix.with_suffix(".json")
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    sidecar.write_text(
        json.dumps({"width_scale": width_scale, "cases": sidecar_cases}, indent=2),
        encoding="utf-8",
    )
    return {"png": str(png), "svg": str(svg), "data_json": str(sidecar), "cases": sidecar_cases}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot generalized-chain snapshots for explicit P2U network cases.")
    parser.add_argument("--r-value", type=int, default=64)
    parser.add_argument("--label", default=None)
    parser.add_argument("--final-gpkg", type=Path, default=None)
    parser.add_argument("--metadata-json", type=Path, default=None)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=OUTPUT_DIR / "p2u_original_vs_source_repaired_R64_generalized_chains",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    final_gpkg = args.final_gpkg or case_paths(args.r_value, "exact")["final_gpkg"]
    metadata = args.metadata_json or case_paths(args.r_value, "exact")["final_metadata"]
    label = args.label or f"Source-repaired R{args.r_value}"
    cases = [
        _original_snapshot_tables(),
        _final_snapshot_tables(final_gpkg, metadata, f"R{args.r_value}", label),
    ]
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    print(json.dumps(plot_cases(cases, output_prefix=args.output_prefix), indent=2))


if __name__ == "__main__":
    main()
