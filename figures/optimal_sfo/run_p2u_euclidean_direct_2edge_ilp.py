from __future__ import annotations

import argparse
import html
import json
import sys
import time
from pathlib import Path
from typing import Any

import geopandas as gpd
import gurobipy as gp
import networkx as nx
import numpy as np
from gurobipy import GRB
from scipy.spatial import Delaunay, QhullError, cKDTree
from shapely.geometry import LineString

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.run_p2u_euclidean_equal_chains import (  # noqa: E402
    INPUT_GPKG,
    load_p2u_terminal_points,
)


OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo" / "euclidean_direct_ilp"
OUTPUT_GPKG = OUTPUT_DIR / "p2u_euclidean_direct_2edge_ilp_3857.gpkg"
OUTPUT_QGS = OUTPUT_DIR / "p2u_euclidean_direct_2edge_ilp.qgs"
SUMMARY_JSON = OUTPUT_DIR / "p2u_euclidean_direct_2edge_ilp_summary.json"
SUMMARY_MD = OUTPUT_DIR / "p2u_euclidean_direct_2edge_ilp_summary.md"


def edge_key(u: int, v: int) -> tuple[int, int]:
    return tuple(sorted((int(u), int(v))))


def cycle_rank(graph: nx.Graph) -> int:
    return graph.number_of_edges() - graph.number_of_nodes() + nx.number_connected_components(graph)


def build_candidate_graph(
    terminals: gpd.GeoDataFrame,
    *,
    k_nearest: int,
    include_delaunay: bool,
    warm_start_gpkg: Path | None = None,
) -> nx.Graph:
    points = np.column_stack([terminals.geometry.x.to_numpy(), terminals.geometry.y.to_numpy()]).astype(float)
    graph = nx.Graph()
    for idx, row in terminals.iterrows():
        graph.add_node(
            int(idx),
            point_index=int(row.point_index),
            terminal_id=str(row.terminal_id),
            kind=str(row.kind),
            source_count=int(row.source_count),
            transformer_count=int(row.transformer_count),
            size_kva=float(row.size_kva),
            nominal_voltage_kv=float(row.nominal_voltage_kv),
            pos=points[int(idx)],
        )

    edges: set[tuple[int, int]] = set()
    if k_nearest > 0:
        tree = cKDTree(points)
        _, neighbors = tree.query(points, k=min(k_nearest + 1, len(points)))
        for i, row in enumerate(np.atleast_2d(neighbors)):
            for j in row:
                j = int(j)
                if i != j:
                    edges.add(edge_key(i, j))

    delaunay_status = "skipped"
    if include_delaunay:
        try:
            tri = Delaunay(points)
            for simplex in tri.simplices:
                for i in range(len(simplex)):
                    for j in range(i + 1, len(simplex)):
                        edges.add(edge_key(int(simplex[i]), int(simplex[j])))
            delaunay_status = "ok"
        except QhullError as exc:
            delaunay_status = f"qhull_error: {exc.__class__.__name__}"

    for u, v in sorted(edges):
        length = float(np.linalg.norm(points[u] - points[v]))
        graph.add_edge(u, v, length_m=length, length=length)

    warm_start_edges = load_warm_start_edges(warm_start_gpkg, graph) if warm_start_gpkg is not None else set()
    for u, v in warm_start_edges:
        if not graph.has_edge(u, v):
            length = float(np.linalg.norm(points[u] - points[v]))
            graph.add_edge(u, v, length_m=length, length=length, warm_start_only=True)

    graph.graph["candidate_summary"] = {
        "k_nearest": int(k_nearest),
        "include_delaunay": bool(include_delaunay),
        "delaunay_status": delaunay_status,
        "warm_start_gpkg": str(warm_start_gpkg) if warm_start_gpkg is not None else None,
        "warm_start_edges": int(len(warm_start_edges)),
        "candidate_edges": int(graph.number_of_edges()),
        "candidate_components": int(nx.number_connected_components(graph)),
        "candidate_bridges": int(len(list(nx.bridges(graph))) if nx.is_connected(graph) else -1),
    }
    graph.graph["warm_start_edges"] = warm_start_edges
    return graph


def load_warm_start_edges(gpkg: Path | None, graph: nx.Graph) -> set[tuple[int, int]]:
    if gpkg is None:
        return set()
    terminal_to_node = {str(data["terminal_id"]): int(node) for node, data in graph.nodes(data=True)}
    layer_candidates = [
        "euclidean_edges",
        "direct_ilp_edges",
        "solution_edges",
        "final_edges",
    ]
    last_error: Exception | None = None
    edges = None
    for layer in layer_candidates:
        try:
            edges = gpd.read_file(gpkg, layer=layer)
            break
        except Exception as exc:  # noqa: BLE001 - try known layer names and fail below with context
            last_error = exc
    if edges is None:
        raise ValueError(f"could not read a recognized edge layer from {gpkg}: {last_error}")
    result: set[tuple[int, int]] = set()
    for _, row in edges.iterrows():
        if "terminal_a" in row and "terminal_b" in row:
            a = str(row.terminal_a)
            b = str(row.terminal_b)
        elif "original_terminal_a" in row and "original_terminal_b" in row:
            a = str(row.original_terminal_a)
            b = str(row.original_terminal_b)
        else:
            continue
        if a not in terminal_to_node or b not in terminal_to_node:
            continue
        u = terminal_to_node[a]
        v = terminal_to_node[b]
        if u != v:
            result.add(edge_key(u, v))
    if not result:
        raise ValueError(f"warm-start file {gpkg} did not contain any usable terminal edge")
    return result


def solve_direct_min_2edge(
    graph: nx.Graph,
    *,
    max_redundancy: int | None,
    time_limit: float,
    mip_gap: float,
    threads: int,
    cut_mode: str,
    max_cut_rounds: int,
) -> tuple[nx.Graph | None, dict[str, Any]]:
    if cut_mode not in {"callback", "iterative"}:
        raise ValueError("cut_mode must be 'callback' or 'iterative'")
    if not nx.is_connected(graph):
        raise ValueError("candidate graph is not connected; increase --k-nearest")
    candidate_bridges = list(nx.bridges(graph))
    if candidate_bridges:
        raise ValueError(f"candidate graph has {len(candidate_bridges)} bridges; increase --k-nearest")

    nodes = sorted(graph.nodes)
    edges = sorted(edge_key(u, v) for u, v in graph.edges)
    incident = {node: [] for node in nodes}
    for u, v in edges:
        incident[u].append((u, v))
        incident[v].append((u, v))

    model = gp.Model("p2u_euclidean_direct_2edge")
    model.Params.OutputFlag = 1
    model.Params.MIPFocus = 1
    model.Params.MIPGap = mip_gap
    model.Params.Threads = threads
    model.Params.Presolve = 2
    model.Params.Heuristics = 0.2
    if time_limit > 0:
        model.Params.TimeLimit = time_limit
    if cut_mode == "callback":
        model.Params.LazyConstraints = 1

    x = model.addVars(edges, vtype=GRB.BINARY, name="x")
    warm_start_edges = set(graph.graph.get("warm_start_edges", set()))
    for edge in edges:
        x[edge].Start = 1.0 if edge in warm_start_edges else 0.0
    for node in nodes:
        model.addConstr(gp.quicksum(x[e] for e in incident[node]) >= 2, name=f"deg2[{node}]")

    if max_redundancy is not None:
        model.addConstr(
            gp.quicksum(x[e] for e in edges) <= len(nodes) - 1 + int(max_redundancy),
            name="max_edge_budget",
        )

    model.setObjective(gp.quicksum(float(graph.edges[e]["length_m"]) * x[e] for e in edges), GRB.MINIMIZE)

    start = time.time()
    cut_rounds = 0
    added_cuts = 0
    selected_edges: list[tuple[int, int]] = []
    stop_reason = "not_started"

    def violated_cutsets(selected: list[tuple[int, int]]) -> list[list[tuple[int, int]]]:
        solution = nx.Graph()
        solution.add_nodes_from(nodes)
        solution.add_edges_from(selected)
        cutsets: list[list[tuple[int, int]]] = []
        comps = list(nx.connected_components(solution))
        if len(comps) > 1:
            for comp in comps:
                crossing = [e for e in edges if (e[0] in comp) ^ (e[1] in comp)]
                if crossing:
                    cutsets.append(crossing)
            return cutsets
        for bridge in nx.bridges(solution):
            bridge_key = edge_key(*bridge)
            tmp = solution.copy()
            tmp.remove_edge(*bridge_key)
            for comp in nx.connected_components(tmp):
                crossing = [e for e in edges if (e[0] in comp) ^ (e[1] in comp)]
                if crossing:
                    cutsets.append(crossing)
        return cutsets

    if cut_mode == "callback":
        callback_stats = {"rounds": 0, "cuts": 0}

        def callback(cb_model, where):
            if where != GRB.Callback.MIPSOL:
                return
            vals = cb_model.cbGetSolution(x)
            incumbent = [e for e in edges if vals[e] > 0.5]
            cutsets = violated_cutsets(incumbent)
            if not cutsets:
                return
            callback_stats["rounds"] += 1
            callback_stats["cuts"] += len(cutsets)
            for crossing in cutsets:
                cb_model.cbLazy(gp.quicksum(x[e] for e in crossing) >= 2)

        model.optimize(callback)
        cut_rounds = int(callback_stats["rounds"])
        added_cuts = int(callback_stats["cuts"])
        if model.SolCount:
            selected_edges = [e for e in edges if x[e].X > 0.5]
            stop_reason = "2edge_connected" if not violated_cutsets(selected_edges) else "violated_cuts_remaining"
        else:
            stop_reason = "no_solution"
    else:
        stop_reason = "max_cut_rounds"
        for round_idx in range(1, max_cut_rounds + 1):
            if time_limit > 0:
                remaining = time_limit - (time.time() - start)
                if remaining <= 0:
                    stop_reason = "global_time_limit"
                    break
                model.Params.TimeLimit = remaining
            model.optimize()
            cut_rounds = round_idx
            if not model.SolCount:
                stop_reason = "no_solution"
                break
            selected_edges = [e for e in edges if x[e].X > 0.5]
            cutsets = violated_cutsets(selected_edges)
            if not cutsets:
                stop_reason = "2edge_connected"
                break
            for crossing in cutsets:
                model.addConstr(gp.quicksum(x[e] for e in crossing) >= 2)
            added_cuts += len(cutsets)
            model.update()

    elapsed = time.time() - start
    summary = {
        "status": int(model.Status),
        "status_name": gurobi_status_name(model.Status),
        "stop_reason": stop_reason,
        "runtime_s": float(elapsed),
        "cut_mode": cut_mode,
        "cut_rounds": int(cut_rounds),
        "added_cuts": int(added_cuts),
        "input_nodes": int(graph.number_of_nodes()),
        "input_edges": int(graph.number_of_edges()),
        "max_redundancy_constraint": max_redundancy,
        "time_limit_s": float(time_limit),
        "threads": int(threads),
        "requested_mip_gap": float(mip_gap),
        "warm_start_edges_used": int(len(warm_start_edges)),
        **graph.graph.get("candidate_summary", {}),
    }
    if not model.SolCount:
        return None, summary

    solution = nx.Graph()
    solution.add_nodes_from((node, graph.nodes[node]) for node in nodes)
    for u, v in selected_edges:
        solution.add_edge(u, v, **graph.edges[u, v])
    bridges = list(nx.bridges(solution)) if nx.is_connected(solution) else []
    summary.update(
        {
            "objective_length_m": float(model.ObjVal),
            "best_bound": float(model.ObjBound),
            "mip_gap": float(model.MIPGap),
            "solution_nodes": int(solution.number_of_nodes()),
            "solution_edges": int(solution.number_of_edges()),
            "solution_cycle_rank": int(cycle_rank(solution)),
            "solution_connected": bool(nx.is_connected(solution)),
            "solution_bridge_count": int(len(bridges)),
            "solution_is_2edge_connected": bool(nx.is_connected(solution) and len(bridges) == 0),
        }
    )
    return solution, summary


def gurobi_status_name(status: int) -> str:
    names = {
        GRB.LOADED: "LOADED",
        GRB.OPTIMAL: "OPTIMAL",
        GRB.INFEASIBLE: "INFEASIBLE",
        GRB.TIME_LIMIT: "TIME_LIMIT",
        GRB.INTERRUPTED: "INTERRUPTED",
        GRB.SUBOPTIMAL: "SUBOPTIMAL",
        GRB.INF_OR_UNBD: "INF_OR_UNBD",
        GRB.UNBOUNDED: "UNBOUNDED",
    }
    return names.get(status, f"STATUS_{status}")


def mark_normally_open_ties(graph: nx.Graph) -> None:
    tree_edges: set[tuple[int, int]] = set()
    for comp in nx.connected_components(graph):
        subgraph = graph.subgraph(comp)
        for u, v, _ in nx.minimum_spanning_edges(subgraph, data=True, weight="length_m"):
            tree_edges.add(edge_key(u, v))
    for u, v, data in graph.edges(data=True):
        is_tie = edge_key(u, v) not in tree_edges
        data["is_tie"] = bool(is_tie)
        data["is_switch"] = bool(is_tie)
        data["normally_closed"] = not is_tie


def graph_to_layers(graph: nx.Graph, terminals: gpd.GeoDataFrame) -> dict[str, gpd.GeoDataFrame]:
    crs = terminals.crs
    mark_normally_open_ties(graph)
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
        pos_u = np.asarray(graph.nodes[u]["pos"], dtype=float)
        pos_v = np.asarray(graph.nodes[v]["pos"], dtype=float)
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
                "geometry": LineString([tuple(pos_u), tuple(pos_v)]),
            }
        )
    edges = gpd.GeoDataFrame(edge_rows, geometry="geometry", crs=crs)
    return {
        "direct_ilp_nodes": nodes,
        "direct_ilp_transformers": nodes[nodes["kind"] == "transformer"].copy(),
        "direct_ilp_sources": nodes[nodes["kind"] == "source"].copy(),
        "direct_ilp_edges": edges,
        "direct_ilp_normally_closed_edges": edges[edges["normally_closed"]].copy(),
        "direct_ilp_normally_open_ties": edges[~edges["normally_closed"]].copy(),
    }


def write_outputs(
    solution: nx.Graph | None,
    terminals: gpd.GeoDataFrame,
    summary: dict[str, Any],
    *,
    output_gpkg: Path,
    output_qgs: Path,
    summary_json: Path,
    summary_md: Path,
) -> None:
    output_gpkg.parent.mkdir(parents=True, exist_ok=True)
    if output_gpkg.exists():
        output_gpkg.unlink()
    if solution is not None:
        for layer, gdf in graph_to_layers(solution, terminals).items():
            gdf.to_file(output_gpkg, layer=layer, driver="GPKG")
        write_qgis_project(output_gpkg, output_qgs)
        summary["output_gpkg"] = str(output_gpkg)
        summary["output_qgs"] = str(output_qgs)
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_md.write_text(summary_markdown(summary), encoding="utf-8")


def summary_markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# P2U Euclidean Direct 2-Edge ILP",
        "",
        "This is a direct sparse-candidate ILP baseline. It is intended to test whether the high cost of the fully 2-connected Euclidean designs comes from 2-connectivity itself or from the equal-chain construction.",
        "",
        f"- Status: `{summary.get('status_name')}`",
        f"- Stop reason: `{summary.get('stop_reason')}`",
        f"- Runtime: `{summary.get('runtime_s', 0):.2f}` s",
        f"- Candidate k-nearest: `{summary.get('k_nearest')}`",
        f"- Delaunay: `{summary.get('delaunay_status')}`",
        f"- Warm start edges: `{summary.get('warm_start_edges')}`",
        f"- Candidate nodes/edges: `{summary.get('input_nodes')}` / `{summary.get('input_edges')}`",
        f"- Candidate bridges: `{summary.get('candidate_bridges')}`",
        f"- Solution nodes/edges: `{summary.get('solution_nodes')}` / `{summary.get('solution_edges')}`",
        f"- Solution physical R: `{summary.get('solution_cycle_rank')}`",
        f"- Solution connected: `{summary.get('solution_connected')}`",
        f"- Solution bridge count: `{summary.get('solution_bridge_count')}`",
        f"- Solution 2-edge-connected: `{summary.get('solution_is_2edge_connected')}`",
        f"- Objective length: `{summary.get('objective_length_m', 0):.3f}` m",
        f"- Best bound: `{summary.get('best_bound', 0):.3f}`",
        f"- MIP gap: `{summary.get('mip_gap')}`",
        f"- Max redundancy constraint: `{summary.get('max_redundancy_constraint')}`",
        f"- Cut mode: `{summary.get('cut_mode')}`",
        f"- Cut rounds: `{summary.get('cut_rounds')}`",
        f"- Added cuts: `{summary.get('added_cuts')}`",
        f"- Threads: `{summary.get('threads')}`",
        f"- Requested MIPGap: `{summary.get('requested_mip_gap')}`",
    ]
    return "\n".join(lines) + "\n"


def write_qgis_project(gpkg_path: Path, qgs_path: Path) -> None:
    specs = [
        ("direct_ilp_normally_closed_edges", "Direct ILP normally closed edges", "direct_closed_1f0ed05f", "Line", line_symbol("142,52,142,255", "0.45")),
        ("direct_ilp_normally_open_ties", "Direct ILP normally open ties", "direct_ties_c539017d", "Line", line_symbol("245,157,35,255", "0.85", "dash")),
        ("direct_ilp_transformers", "Transformer points", "direct_transformers_b1648e42", "Point", point_symbol("255,255,255,0", "55,55,55,255", "0.75", False)),
        ("direct_ilp_sources", "Source points", "direct_sources_2cbd31ac", "Point", point_symbol("18,18,18,255", "255,255,255,255", "2.4", True)),
    ]
    tree = "\n".join(
        f'      <layer-tree-layer checked="1" id="{layer_id}" name="{html.escape(name)}" source="./{gpkg_path.name}|layername={layer}"/>'
        for layer, name, layer_id, _, _ in specs
    )
    layers = "\n".join(
        map_layer(gpkg_path.name, layer, name, layer_id, geometry, renderer)
        for layer, name, layer_id, geometry, renderer in specs
    )
    title = gpkg_path.stem.replace("_3857", "").replace("_", " ")
    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="{html.escape(title)}">
  <homePath path="."/>
  <title>{html.escape(title)}</title>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="{html.escape(title)}">
{tree}
    </layer-tree-group>
  </layer-tree-group>
  <projectlayers>
{layers}
  </projectlayers>
</qgis>
"""
    qgs_path.write_text(qgs, encoding="utf-8")


def line_symbol(color: str, width: str, style: str = "solid") -> str:
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


def point_symbol(color: str, outline: str, size: str, fill: bool) -> str:
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


def map_layer(gpkg_name: str, layer: str, name: str, layer_id: str, geometry: str, renderer: str) -> str:
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run direct Euclidean minimum 2-edge-connected ILP on all P2U terminals.")
    parser.add_argument("--input-gpkg", type=Path, default=INPUT_GPKG)
    parser.add_argument("--k-nearest", type=int, default=8)
    parser.add_argument("--no-delaunay", action="store_true")
    parser.add_argument("--warm-start-gpkg", type=Path, default=None)
    parser.add_argument("--max-redundancy", type=int, default=None)
    parser.add_argument("--time-limit", type=float, default=900.0)
    parser.add_argument("--mip-gap", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--cut-mode", choices=["callback", "iterative"], default="callback")
    parser.add_argument("--max-cut-rounds", type=int, default=100)
    parser.add_argument("--output-gpkg", type=Path, default=OUTPUT_GPKG)
    parser.add_argument("--output-qgs", type=Path, default=OUTPUT_QGS)
    parser.add_argument("--summary-json", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--summary-md", type=Path, default=SUMMARY_MD)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    terminals = load_p2u_terminal_points(args.input_gpkg)
    candidate = build_candidate_graph(
        terminals,
        k_nearest=args.k_nearest,
        include_delaunay=not args.no_delaunay,
        warm_start_gpkg=args.warm_start_gpkg,
    )
    solution, summary = solve_direct_min_2edge(
        candidate,
        max_redundancy=args.max_redundancy,
        time_limit=args.time_limit,
        mip_gap=args.mip_gap,
        threads=args.threads,
        cut_mode=args.cut_mode,
        max_cut_rounds=args.max_cut_rounds,
    )
    summary.update(
        {
            "input_gpkg": str(args.input_gpkg),
            "output_gpkg": str(args.output_gpkg),
            "output_qgs": str(args.output_qgs),
        }
    )
    write_outputs(
        solution,
        terminals,
        summary,
        output_gpkg=args.output_gpkg,
        output_qgs=args.output_qgs,
        summary_json=args.summary_json,
        summary_md=args.summary_md,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
