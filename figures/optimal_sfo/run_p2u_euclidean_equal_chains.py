from __future__ import annotations

import argparse
import html
import json
import math
import sys
import time
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
from shapely.geometry import LineString

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from optimal_network import add_chains_to_strc, balanced_kmeans_gurobi, get_optimal_strc_trips, optimize_rel_weight_ratio


INPUT_GPKG = ROOT / "outputs" / "qgis" / "P2U" / "P2U_MV_physical_network_3857.gpkg"
OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo" / "euclidean_equal_chains"
OUTPUT_GPKG = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_3857.gpkg"
OUTPUT_QGS = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50.qgs"
SUMMARY_JSON = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_summary.json"
SUMMARY_MD = OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_summary.md"


def load_p2u_terminal_points(input_gpkg: Path = INPUT_GPKG) -> gpd.GeoDataFrame:
    transformers = gpd.read_file(input_gpkg, layer="mv_transformer_load_points")
    sources = gpd.read_file(input_gpkg, layer="mv_sources_hvmv_substations")
    if transformers.crs != sources.crs:
        sources = sources.to_crs(transformers.crs)
    if transformers.crs is None:
        raise ValueError("input point layers must have a CRS")

    transformer_rows = transformers.copy()
    transformer_rows["terminal_id"] = "T:" + transformer_rows["Node"].astype(str)
    transformer_rows["kind"] = "transformer"
    transformer_rows["source_count"] = 0
    transformer_rows["transformer_count"] = 1
    transformer_rows["size_kva"] = transformer_rows["Size_KVA" if "Size_KVA" in transformer_rows.columns else "Size_kVA" if "Size_kVA" in transformer_rows.columns else "Size_KVA"].astype(float)
    transformer_rows["nominal_voltage_kv"] = transformer_rows["NomV_kV"].astype(float)
    transformer_rows = transformer_rows[
        ["terminal_id", "Node", "kind", "source_count", "transformer_count", "size_kva", "nominal_voltage_kv", "geometry"]
    ]

    source_rows = sources.copy()
    source_rows["terminal_id"] = "S:" + source_rows["Node"].astype(str)
    source_rows["kind"] = "source"
    source_rows["source_count"] = 1
    source_rows["transformer_count"] = 0
    source_rows["size_kva"] = 0.0
    source_rows["nominal_voltage_kv"] = 12.47
    source_rows = source_rows[
        ["terminal_id", "Node", "kind", "source_count", "transformer_count", "size_kva", "nominal_voltage_kv", "geometry"]
    ]

    terminals = gpd.GeoDataFrame(
        pd.concat([transformer_rows, source_rows], ignore_index=True),
        geometry="geometry",
        crs=transformers.crs,
    )
    terminals = terminals.drop_duplicates(subset=["terminal_id"], keep="first").reset_index(drop=True)
    terminals["point_index"] = terminals.index.astype(int)
    return terminals


def run_equal_chain_algorithm(
    terminals: gpd.GeoDataFrame,
    *,
    redundancy: int,
    seed: int,
    kmeans_max_iter: int,
    strc_n_init_iters: int,
    strc_exact_vertices: bool,
    strc_trip_nearest_vertices: int | None,
    chain_n_init_iters: int,
    local_fix_max_changes: int,
    local_fix_max_risk_gain: float,
    debug: bool,
    checkpoint_dir: Path | None,
) -> nx.Graph:
    points = np.column_stack([terminals.geometry.x.to_numpy(), terminals.geometry.y.to_numpy()]).astype(float)
    if redundancy > 1:
        n_structure_nodes = 2 * (redundancy - 1)
        n_chains = 3 * (redundancy - 1)
    elif redundancy == 1:
        n_structure_nodes = 1
        n_chains = 1
    else:
        raise ValueError("redundancy must be positive")

    stage_t0 = time.perf_counter()
    print(
        f"[stage] balanced_kmeans_gurobi: points={len(points)}, chains={n_chains}, max_iter={kmeans_max_iter}",
        flush=True,
    )
    chains, centers = balanced_kmeans_gurobi(
        X=points,
        k=n_chains,
        max_iter=kmeans_max_iter,
        chain_len_sigma=0.0,
        random_state=seed,
    )
    print(f"[stage] balanced_kmeans_gurobi done in {time.perf_counter() - stage_t0:.2f}s", flush=True)
    if checkpoint_dir is not None:
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            checkpoint_dir / f"p2u_euclidean_equal_chains_Rmax{redundancy}_clusters.npz",
            chains=chains,
            centers=centers,
        )

    stage_t0 = time.perf_counter()
    print(
        f"[stage] get_optimal_strc_trips: structure_nodes={n_structure_nodes}, init_iters={strc_n_init_iters}, exact_vertices={strc_exact_vertices}",
        flush=True,
    )
    chosen_trips = get_optimal_strc_trips(
        centers=centers,
        points=points,
        n_nodes=n_structure_nodes,
        strc_n_init_iters=strc_n_init_iters,
        exact_vertices=strc_exact_vertices,
        max_trip_vertices_per_center=strc_trip_nearest_vertices,
        source_node=None,
        debug=debug,
    )
    print(f"[stage] get_optimal_strc_trips done in {time.perf_counter() - stage_t0:.2f}s", flush=True)
    if checkpoint_dir is not None:
        checkpoint_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(chosen_trips, columns=["terminal_a_index", "terminal_b_index", "chain"]).to_csv(
            checkpoint_dir / f"p2u_euclidean_equal_chains_Rmax{redundancy}_structure_trips.csv",
            index=False,
        )

    stage_t0 = time.perf_counter()
    print(
        f"[stage] add_chains_to_strc: trips={len(chosen_trips)}, init_iters={chain_n_init_iters}",
        flush=True,
    )
    graph = add_chains_to_strc(
        points=points,
        chains=chains,
        chosen_trips=chosen_trips,
        max_init_iter=chain_n_init_iters,
        debug=debug,
    )
    print(f"[stage] add_chains_to_strc done in {time.perf_counter() - stage_t0:.2f}s", flush=True)

    terminal_by_index = terminals.set_index("point_index")
    for node in graph.nodes:
        row = terminal_by_index.loc[int(node)]
        graph.nodes[node].update(
            terminal_id=str(row.terminal_id),
            kind=str(row.kind),
            source_count=int(row.source_count),
            transformer_count=int(row.transformer_count),
            size_kva=float(row.size_kva),
            nominal_voltage_kv=float(row.nominal_voltage_kv),
            pos=np.asarray([float(row.geometry.x), float(row.geometry.y)], dtype=float),
        )
    refresh_edge_attributes(graph)

    if local_fix_max_changes > 0:
        before_length = total_graph_length(graph)
        before_risk = internal_chain_risk_proxy(graph)
        before_r = cycle_rank(graph)
        before_bridges = len(list(nx.bridges(graph))) if nx.is_connected(graph) else None
        stage_t0 = time.perf_counter()
        print(
            f"[stage] local fixes: max_changes={local_fix_max_changes}, max_risk_gain={local_fix_max_risk_gain}",
            flush=True,
        )
        graph = optimize_rel_weight_ratio(
            graph,
            max_risk_gain=local_fix_max_risk_gain,
            max_changes=local_fix_max_changes,
            source=None,
            debug=debug,
        )
        refresh_edge_attributes(graph)
        after_r = cycle_rank(graph)
        after_bridges = len(list(nx.bridges(graph))) if nx.is_connected(graph) else None
        if after_r != before_r:
            raise RuntimeError(f"local fixes changed cycle rank from {before_r} to {after_r}")
        if after_bridges not in {0, before_bridges}:
            raise RuntimeError(f"local fixes changed bridge count from {before_bridges} to {after_bridges}")
        graph.graph["local_fix_summary"] = {
            "enabled": True,
            "max_changes": int(local_fix_max_changes),
            "max_risk_gain": float(local_fix_max_risk_gain),
            "runtime_s": float(time.perf_counter() - stage_t0),
            "length_before_m": float(before_length),
            "length_after_m": float(total_graph_length(graph)),
            "length_delta_m": float(total_graph_length(graph) - before_length),
            "length_relative_change": float(total_graph_length(graph) / before_length - 1.0) if before_length else 0.0,
            "internal_chain_risk_before": float(before_risk),
            "internal_chain_risk_after": float(internal_chain_risk_proxy(graph)),
            "cycle_rank_before": int(before_r),
            "cycle_rank_after": int(after_r),
            "bridges_before": before_bridges,
            "bridges_after": after_bridges,
        }
        print(f"[stage] local fixes done in {graph.graph['local_fix_summary']['runtime_s']:.2f}s", flush=True)
    else:
        graph.graph["local_fix_summary"] = {"enabled": False}
    return graph


def refresh_edge_attributes(graph: nx.Graph) -> None:
    for u, v, data in graph.edges(data=True):
        pos_u = np.asarray(graph.nodes[u]["pos"], dtype=float)
        pos_v = np.asarray(graph.nodes[v]["pos"], dtype=float)
        length = float(np.linalg.norm(pos_u - pos_v))
        data["length_m"] = length
        data["length"] = length
        data.setdefault("is_tie", False)
        data.setdefault("is_switch", False)
        data.setdefault("normally_closed", True)
    mark_normally_open_ties(graph)


def mark_normally_open_ties(graph: nx.Graph) -> None:
    tree_edges: set[tuple[int, int]] = set()
    for component in nx.connected_components(graph):
        subgraph = graph.subgraph(component)
        for u, v, _ in nx.minimum_spanning_edges(subgraph, data=True, weight="length_m"):
            tree_edges.add(edge_key(u, v))
    for u, v, data in graph.edges(data=True):
        is_tie = edge_key(u, v) not in tree_edges
        data["is_tie"] = bool(is_tie)
        data["is_switch"] = bool(is_tie)
        data["normally_closed"] = not is_tie


def cycle_rank(graph: nx.Graph) -> int:
    return graph.number_of_edges() - graph.number_of_nodes() + nx.number_connected_components(graph)


def total_graph_length(graph: nx.Graph) -> float:
    return float(sum(float(data.get("length_m", 0.0)) for _, _, data in graph.edges(data=True)))


def internal_chain_risk_proxy(graph: nx.Graph) -> float:
    from optimal_network.reduce_weight import _total_inter_risk

    return float(_total_inter_risk(graph))


def edge_key(u: int, v: int) -> tuple[int, int]:
    return tuple(sorted((int(u), int(v))))


def graph_to_geodataframes(graph: nx.Graph, terminals: gpd.GeoDataFrame) -> dict[str, gpd.GeoDataFrame]:
    crs = terminals.crs
    node_rows = []
    for node, data in graph.nodes(data=True):
        point = terminals.loc[int(node)].geometry
        node_rows.append(
            {
                "point_index": int(node),
                "terminal_id": data["terminal_id"],
                "kind": data["kind"],
                "source_count": int(data["source_count"]),
                "transformer_count": int(data["transformer_count"]),
                "size_kva": float(data["size_kva"]),
                "nominal_voltage_kv": float(data["nominal_voltage_kv"]),
                "geometry": point,
            }
        )
    nodes = gpd.GeoDataFrame(node_rows, geometry="geometry", crs=crs)

    edge_rows = []
    for edge_id, (u, v, data) in enumerate(graph.edges(data=True)):
        pos_u = graph.nodes[u]["pos"]
        pos_v = graph.nodes[v]["pos"]
        edge_rows.append(
            {
                "edge_id": int(edge_id),
                "u_index": int(u),
                "v_index": int(v),
                "terminal_a": graph.nodes[u]["terminal_id"],
                "terminal_b": graph.nodes[v]["terminal_id"],
                "kind_a": graph.nodes[u]["kind"],
                "kind_b": graph.nodes[v]["kind"],
                "length_m": float(data["length_m"]),
                "is_tie": bool(data.get("is_tie", False)),
                "is_switch": bool(data.get("is_switch", False)),
                "normally_closed": bool(data.get("normally_closed", True)),
                "geometry": LineString([pos_u, pos_v]),
            }
        )
    edges = gpd.GeoDataFrame(edge_rows, geometry="geometry", crs=crs)
    return {
        "euclidean_nodes": nodes,
        "euclidean_transformers": nodes[nodes["kind"] == "transformer"].copy(),
        "euclidean_sources": nodes[nodes["kind"] == "source"].copy(),
        "euclidean_edges": edges,
        "euclidean_normally_closed_edges": edges[edges["normally_closed"]].copy(),
        "euclidean_normally_open_ties": edges[~edges["normally_closed"]].copy(),
    }


def write_outputs(
    *,
    tables: dict[str, gpd.GeoDataFrame],
    graph: nx.Graph,
    redundancy: int,
    runtime_s: float,
    output_gpkg: Path,
    output_qgs: Path,
    summary_json: Path,
    summary_md: Path,
    args: argparse.Namespace,
) -> dict:
    output_gpkg.parent.mkdir(parents=True, exist_ok=True)
    if output_gpkg.exists():
        output_gpkg.unlink()
    for layer, gdf in tables.items():
        gdf.to_file(output_gpkg, layer=layer, driver="GPKG")
    write_qgis_project(output_gpkg, output_qgs)

    lengths = [float(data["length_m"]) for _, _, data in graph.edges(data=True)]
    tie_count = sum(1 for _, _, data in graph.edges(data=True) if data.get("is_tie", False))
    summary = {
        "input_gpkg": str(args.input_gpkg),
        "output_gpkg": str(output_gpkg),
        "output_qgs": str(output_qgs),
        "redundancy_requested": int(redundancy),
        "seed": int(args.seed),
        "kmeans_max_iter": int(args.kmeans_max_iter),
        "strc_n_init_iters": int(args.strc_n_init_iters),
        "strc_exact_vertices": bool(args.strc_exact_vertices),
        "strc_trip_nearest_vertices": args.strc_trip_nearest_vertices,
        "chain_n_init_iters": int(args.chain_n_init_iters),
        "local_fix": graph.graph.get("local_fix_summary", {"enabled": False}),
        "runtime_s": float(runtime_s),
        "nodes": int(graph.number_of_nodes()),
        "edges": int(graph.number_of_edges()),
        "connected_components": int(nx.number_connected_components(graph)),
        "cycle_rank_R": int(cycle_rank(graph)),
        "tie_edges": int(tie_count),
        "normally_closed_edges": int(graph.number_of_edges() - tie_count),
        "source_nodes": int(sum(1 for _, data in graph.nodes(data=True) if data.get("kind") == "source")),
        "transformer_nodes": int(sum(1 for _, data in graph.nodes(data=True) if data.get("kind") == "transformer")),
        "total_length_m": float(sum(lengths)),
        "mean_edge_length_m": float(np.mean(lengths)) if lengths else 0.0,
        "max_edge_length_m": float(max(lengths)) if lengths else 0.0,
        "is_connected": bool(nx.is_connected(graph)),
    }
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_md.write_text(summary_to_markdown(summary), encoding="utf-8")
    return summary


def write_qgis_project(gpkg_path: Path, qgs_path: Path) -> None:
    specs = [
        ("euclidean_normally_closed_edges", "Euclidean normally closed edges", "closed_edges_4d739d4e", "Line", _line_symbol("203,64,60,255", "0.45")),
        ("euclidean_normally_open_ties", "Euclidean normally open ties", "open_ties_ad74edf3", "Line", _line_symbol("245,157,35,255", "0.85", "dash")),
        ("euclidean_transformers", "Transformer points", "transformers_9cc0a42d", "Point", _point_symbol("255,255,255,0", "55,55,55,255", "0.75", False)),
        ("euclidean_sources", "Source points", "sources_2d3d5024", "Point", _point_symbol("18,18,18,255", "255,255,255,255", "2.4", True)),
    ]
    trees = []
    layers = []
    for layer, name, layer_id, geometry, renderer in specs:
        trees.append(
            f'      <layer-tree-layer checked="1" id="{layer_id}" name="{html.escape(name)}" source="./{gpkg_path.name}|layername={layer}"/>'
        )
        layers.append(_map_layer(gpkg_path.name, layer, name, layer_id, geometry, renderer))

    project_title = gpkg_path.stem.replace("_3857", "").replace("_", " ")
    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="{html.escape(project_title)}">
  <homePath path="."/>
  <title>{html.escape(project_title)}</title>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="{html.escape(project_title)}">
{chr(10).join(trees)}
    </layer-tree-group>
  </layer-tree-group>
  <projectlayers>
{chr(10).join(layers)}
  </projectlayers>
</qgis>
"""
    qgs_path.write_text(qgs, encoding="utf-8")


def _line_symbol(color: str, width: str, style: str = "solid") -> str:
    return f"""      <renderer-v2 type="singleSymbol">
        <symbols>
          <symbol alpha="1" type="line" name="0">
            <layer enabled="1" class="SimpleLine">
              <Option type="Map">
                <Option name="line_color" type="QString" value="{color}"/>
                <Option name="line_style" type="QString" value="{style}"/>
                <Option name="line_width" type="QString" value="{width}"/>
                <Option name="line_width_unit" type="QString" value="MM"/>
                <Option name="capstyle" type="QString" value="round"/>
                <Option name="joinstyle" type="QString" value="round"/>
              </Option>
            </layer>
          </symbol>
        </symbols>
      </renderer-v2>"""


def _point_symbol(color: str, outline: str, size: str, fill: bool) -> str:
    fill_color = color if fill else "255,255,255,0"
    return f"""      <renderer-v2 type="singleSymbol">
        <symbols>
          <symbol alpha="1" type="marker" name="0">
            <layer enabled="1" class="SimpleMarker">
              <Option type="Map">
                <Option name="color" type="QString" value="{fill_color}"/>
                <Option name="outline_color" type="QString" value="{outline}"/>
                <Option name="outline_width" type="QString" value="0.30"/>
                <Option name="outline_width_unit" type="QString" value="MM"/>
                <Option name="size" type="QString" value="{size}"/>
                <Option name="size_unit" type="QString" value="MM"/>
              </Option>
            </layer>
          </symbol>
        </symbols>
      </renderer-v2>"""


def _map_layer(gpkg_name: str, layer: str, name: str, layer_id: str, geometry: str, renderer: str) -> str:
    return f"""  <maplayer type="vector" geometry="{geometry}">
    <id>{layer_id}</id>
    <datasource>./{gpkg_name}|layername={layer}</datasource>
    <layername>{html.escape(name)}</layername>
    <srs>
      <spatialrefsys>
        <authid>EPSG:3857</authid>
        <description>WGS 84 / Pseudo-Mercator</description>
        <projectionacronym>merc</projectionacronym>
        <ellipsoidacronym>EPSG:7030</ellipsoidacronym>
        <geographicflag>false</geographicflag>
      </spatialrefsys>
    </srs>
    <provider encoding="UTF-8">ogr</provider>
{renderer}
  </maplayer>"""


def summary_to_markdown(summary: dict) -> str:
    local_fix = summary.get("local_fix", {"enabled": False})
    description = (
        "This is a first-stage Euclidean network. Edges are straight-line geometric connections; no road embedding has been applied."
        if local_fix.get("enabled")
        else "This is a first-stage Euclidean network. Edges are straight-line geometric connections; no road embedding or cost fine tuning has been applied."
    )
    local_fix_lines = []
    if local_fix.get("enabled"):
        local_fix_lines = [
            "",
            "## Local Fixes",
            "",
            f"- Max changes: `{local_fix['max_changes']}`",
            f"- Max internal-risk gain: `{local_fix['max_risk_gain']}`",
            f"- Length before/after: `{local_fix['length_before_m']:.3f}` m / `{local_fix['length_after_m']:.3f}` m",
            f"- Relative length change: `{local_fix['length_relative_change']:.6g}`",
            f"- Internal risk proxy before/after: `{local_fix['internal_chain_risk_before']:.6g}` / `{local_fix['internal_chain_risk_after']:.6g}`",
        ]
    return "\n".join(
        [
            "# P2U Euclidean Equal-Chain Network",
            "",
            description,
            "",
            f"- Requested R: `{summary['redundancy_requested']}`",
            f"- Realized R: `{summary['cycle_rank_R']}`",
            f"- Nodes/edges: `{summary['nodes']}` / `{summary['edges']}`",
            f"- Sources/transformers: `{summary['source_nodes']}` / `{summary['transformer_nodes']}`",
            f"- Tie edges: `{summary['tie_edges']}`",
            f"- Connected: `{summary['is_connected']}`",
            f"- Total Euclidean length: `{summary['total_length_m']:.3f}` m",
            f"- Mean/max edge length: `{summary['mean_edge_length_m']:.3f}` m / `{summary['max_edge_length_m']:.3f}` m",
            f"- Runtime: `{summary['runtime_s']:.2f}` s",
            "",
            "Normally open ties are defined for display as edges outside a minimum-length spanning tree of the Euclidean graph.",
            *local_fix_lines,
        ]
    ) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Euclidean equal-chain optimal algorithm on all P2U MV terminal points.")
    parser.add_argument("--input-gpkg", type=Path, default=INPUT_GPKG)
    parser.add_argument("--redundancy", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--kmeans-max-iter", type=int, default=7)
    parser.add_argument("--strc-n-init-iters", type=int, default=2)
    parser.add_argument("--strc-exact-vertices", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--strc-trip-nearest-vertices",
        type=int,
        default=None,
        help="Limit each chain center to pairs among its k nearest structure vertices.",
    )
    parser.add_argument("--chain-n-init-iters", type=int, default=2)
    parser.add_argument("--local-fix-max-changes", type=int, default=0)
    parser.add_argument("--local-fix-max-risk-gain", type=float, default=0.1)
    parser.add_argument("--output-gpkg", type=Path, default=OUTPUT_GPKG)
    parser.add_argument("--output-qgs", type=Path, default=OUTPUT_QGS)
    parser.add_argument("--summary-json", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--summary-md", type=Path, default=SUMMARY_MD)
    parser.add_argument("--checkpoint-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    _resolve_default_output_paths(args)
    return args


def _resolve_default_output_paths(args: argparse.Namespace) -> None:
    suffix = f"Rmax{args.redundancy}"
    if args.local_fix_max_changes > 0:
        suffix += "_localfix"
    stem = f"p2u_euclidean_equal_chains_{suffix}"
    defaults = {
        "output_gpkg": OUTPUT_GPKG,
        "output_qgs": OUTPUT_QGS,
        "summary_json": SUMMARY_JSON,
        "summary_md": SUMMARY_MD,
    }
    replacements = {
        "output_gpkg": OUTPUT_DIR / f"{stem}_3857.gpkg",
        "output_qgs": OUTPUT_DIR / f"{stem}.qgs",
        "summary_json": OUTPUT_DIR / f"{stem}_summary.json",
        "summary_md": OUTPUT_DIR / f"{stem}_summary.md",
    }
    for attr, default in defaults.items():
        if getattr(args, attr) == default:
            setattr(args, attr, replacements[attr])


def main() -> None:
    args = parse_args()
    terminals = load_p2u_terminal_points(args.input_gpkg)
    t0 = time.perf_counter()
    graph = run_equal_chain_algorithm(
        terminals,
        redundancy=args.redundancy,
        seed=args.seed,
        kmeans_max_iter=args.kmeans_max_iter,
        strc_n_init_iters=args.strc_n_init_iters,
        strc_exact_vertices=args.strc_exact_vertices,
        strc_trip_nearest_vertices=args.strc_trip_nearest_vertices,
        chain_n_init_iters=args.chain_n_init_iters,
        local_fix_max_changes=args.local_fix_max_changes,
        local_fix_max_risk_gain=args.local_fix_max_risk_gain,
        debug=args.debug,
        checkpoint_dir=args.checkpoint_dir,
    )
    runtime = time.perf_counter() - t0
    tables = graph_to_geodataframes(graph, terminals)
    summary = write_outputs(
        tables=tables,
        graph=graph,
        redundancy=args.redundancy,
        runtime_s=runtime,
        output_gpkg=args.output_gpkg,
        output_qgs=args.output_qgs,
        summary_json=args.summary_json,
        summary_md=args.summary_md,
        args=args,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
