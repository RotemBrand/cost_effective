from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import LineString


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.p2u_full_mv_reliability import (  # noqa: E402
    aggregate_lv_consumers_to_transformers,
    build_full_mv_graph,
)


DATA_DIR = ROOT / "data" / "power" / "better_grids" / "SFO" / "P2U"
OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo"
BACKBONE_SUMMARY = OUTPUT_DIR / "p2u_ilp_2edge_solution_summary.json"
FINAL_NETWORK_GPKG = OUTPUT_DIR / "p2u_final_network_3857.gpkg"
FINAL_NETWORK_METADATA = OUTPUT_DIR / "p2u_final_network_metadata.json"

SOURCE_CRS = "EPSG:32610"
TARGET_CRS = "EPSG:3857"
SUPER_SOURCE = "__SOURCE__"


def clean(value) -> str:
    return str(value).strip()


def edge_key(u: str, v: str) -> tuple[str, str]:
    return tuple(sorted((str(u), str(v))))


def read_p2u_layer(name: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(DATA_DIR / f"{name}.shp")
    if gdf.crs is None:
        gdf = gdf.set_crs(SOURCE_CRS)
    return gdf


def nearest_street_nodes(points: gpd.GeoDataFrame, street_nodes: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    points = points.to_crs(street_nodes.crs)
    point_geom = points.geometry.centroid
    node_geom = street_nodes.geometry
    tree = cKDTree(np.column_stack([node_geom.x.to_numpy(), node_geom.y.to_numpy()]))
    distances, indices = tree.query(np.column_stack([point_geom.x.to_numpy(), point_geom.y.to_numpy()]), k=1)
    return street_nodes.iloc[indices]["Node"].to_numpy(), distances


def build_road_graph(
    branches: gpd.GeoDataFrame,
    street_nodes: gpd.GeoDataFrame,
) -> tuple[nx.Graph, dict[str, tuple[float, float]]]:
    node_pos = {
        clean(row.Node): (float(row.geometry.x), float(row.geometry.y))
        for _, row in street_nodes.iterrows()
    }
    graph = nx.Graph()
    for _, row in branches.iterrows():
        u = clean(row.Node_A)
        v = clean(row.Node_B)
        if u == v:
            continue
        length_m = float(row.geometry.length)
        if graph.has_edge(u, v):
            graph.edges[u, v]["length_m"] = min(graph.edges[u, v]["length_m"], length_m)
        else:
            graph.add_edge(u, v, length_m=length_m)
    return graph, node_pos


def line_from_road_nodes(nodes: list[str], node_pos: dict[str, tuple[float, float]]) -> LineString:
    coords = []
    for node in nodes:
        coord = node_pos.get(node)
        if coord is None:
            continue
        if not coords or coords[-1] != coord:
            coords.append(coord)
    if len(coords) == 1:
        coords.append(coords[0])
    return LineString(coords)


def terminal_road_node(terminal: str) -> str:
    if ":" not in terminal:
        return terminal
    return terminal.split(":", 1)[1]


def sequence_road_nodes(sequence: str) -> list[str]:
    return [terminal_road_node(token) for token in str(sequence).split("|") if token]


def source_contracted_graph(graph: nx.Graph, sources: set[str]) -> nx.Graph:
    contracted = nx.Graph()
    for node, data in graph.nodes(data=True):
        target = SUPER_SOURCE if node in sources else node
        if target not in contracted:
            contracted.add_node(target, **data)
        else:
            contracted.nodes[target]["weight"] = float(contracted.nodes[target].get("weight", 0.0)) + float(
                data.get("weight", 0.0)
            )
    for u, v, data in graph.edges(data=True):
        cu = SUPER_SOURCE if u in sources else u
        cv = SUPER_SOURCE if v in sources else v
        if cu == cv:
            continue
        length_m = float(data.get("length_m", 1.0))
        if contracted.has_edge(cu, cv) and contracted.edges[cu, cv]["length_m"] <= length_m:
            continue
        contracted.add_edge(cu, cv, **data)
    return contracted


def cycle_rank(graph: nx.Graph) -> int:
    return graph.number_of_edges() - graph.number_of_nodes() + nx.number_connected_components(graph)


def percent(part: float, whole: float) -> float:
    return 0.0 if whole == 0 else 100.0 * float(part) / float(whole)


def load_transformer_demand_by_node() -> tuple[dict[str, dict[str, float]], dict]:
    mv_graph, _, _ = build_full_mv_graph()
    lv_assignment = aggregate_lv_consumers_to_transformers(mv_graph)
    loads = {}
    for node, data in mv_graph.nodes(data=True):
        loads[clean(node)] = {
            "demand_kw": float(data.get("demand_kw", 0.0)),
            "consumer_count": float(data.get("consumer_count", 0.0)),
            "yearly_kwh": float(data.get("yearly_kwh", 0.0)),
        }
    return loads, lv_assignment


def transformers_by_street_node(street_nodes: gpd.GeoDataFrame) -> tuple[dict[str, dict], gpd.GeoDataFrame, dict]:
    transformers = read_p2u_layer("DistribTransf_N").copy()
    transformer_loads, lv_assignment = load_transformer_demand_by_node()
    street_nodes_for_join = street_nodes.copy()
    street_nodes_for_join["Node"] = street_nodes_for_join["Node"].map(clean)
    nearest, distances = nearest_street_nodes(transformers, street_nodes_for_join)

    transformers["transformer"] = transformers["Node"].map(clean)
    transformers["street_node"] = list(map(clean, nearest))
    transformers["nearest_street_distance_m"] = distances.astype(float)
    transformers["size_kva"] = transformers["Size_kVA"].astype(float)
    transformers["demand_kw"] = transformers["transformer"].map(lambda n: transformer_loads.get(n, {}).get("demand_kw", 0.0))
    transformers["consumer_count"] = transformers["transformer"].map(
        lambda n: transformer_loads.get(n, {}).get("consumer_count", 0.0)
    )
    transformers["yearly_kwh"] = transformers["transformer"].map(
        lambda n: transformer_loads.get(n, {}).get("yearly_kwh", 0.0)
    )

    grouped: dict[str, dict] = {}
    for street_node, group in transformers.groupby("street_node"):
        grouped[clean(street_node)] = {
            "transformer_count": int(len(group)),
            "size_kva": float(group["size_kva"].sum()),
            "demand_kw": float(group["demand_kw"].sum()),
            "consumer_count": float(group["consumer_count"].sum()),
            "yearly_kwh": float(group["yearly_kwh"].sum()),
            "transformers": "|".join(group["transformer"].astype(str).sort_values().tolist()),
        }
    return grouped, transformers, lv_assignment


def load_backbone_solution(
    summary_path: Path = BACKBONE_SUMMARY,
) -> tuple[Path, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    gpkg = Path(summary["output_gpkg"])
    return (
        gpkg,
        gpd.read_file(gpkg, layer="solution_edges"),
        gpd.read_file(gpkg, layer="solution_transformer_nodes"),
        gpd.read_file(gpkg, layer="solution_source_nodes"),
    )


def solution_node_rows(
    solution_nodes: gpd.GeoDataFrame,
    terminal_class: str,
    transformer_by_road: dict[str, dict],
    node_pos: dict[str, tuple[float, float]],
) -> list[dict]:
    rows = []
    for _, row in solution_nodes.iterrows():
        road_node = clean(row.road_node)
        transformer_data = transformer_by_road.get(road_node, {})
        rows.append(
            {
                "terminal_id": str(row.terminal_id),
                "road_node": road_node,
                "kind": str(row.kind),
                "terminal_class": terminal_class,
                "transformer_count": int(transformer_data.get("transformer_count", row.get("transformer_count", 0))),
                "source_count": int(row.get("source_count", 0)),
                "size_kva": float(transformer_data.get("size_kva", row.get("size_kva", 0.0))),
                "demand_kw": float(transformer_data.get("demand_kw", 0.0)),
                "consumer_count": float(transformer_data.get("consumer_count", 0.0)),
                "yearly_kwh": float(transformer_data.get("yearly_kwh", 0.0)),
                "geometry": line_from_road_nodes([road_node], node_pos).centroid,
            }
        )
    return rows


def edge_internal_load(sequence: str, transformer_by_road: dict[str, dict]) -> dict[str, float | int]:
    internal_nodes = sequence_road_nodes(sequence)[1:-1]
    return {
        "edge_transformer_count": int(sum(transformer_by_road.get(node, {}).get("transformer_count", 0) for node in internal_nodes)),
        "edge_size_kva": float(sum(transformer_by_road.get(node, {}).get("size_kva", 0.0) for node in internal_nodes)),
        "edge_demand_kw": float(sum(transformer_by_road.get(node, {}).get("demand_kw", 0.0) for node in internal_nodes)),
        "edge_consumer_count": float(sum(transformer_by_road.get(node, {}).get("consumer_count", 0.0) for node in internal_nodes)),
        "edge_yearly_kwh": float(sum(transformer_by_road.get(node, {}).get("yearly_kwh", 0.0) for node in internal_nodes)),
    }


def graph_from_final_tables(nodes: gpd.GeoDataFrame, edges: gpd.GeoDataFrame) -> nx.Graph:
    graph = nx.Graph()
    for _, row in nodes.iterrows():
        graph.add_node(
            str(row.terminal_id),
            kind=str(row.kind),
            terminal_class=str(row.terminal_class),
            road_node=str(row.road_node),
            transformer_count=int(row.transformer_count),
            source_count=int(row.source_count),
            size_kva=float(row.size_kva),
            weight=float(row.demand_kw),
            demand_kw=float(row.demand_kw),
            consumer_count=float(row.consumer_count),
            yearly_kwh=float(row.yearly_kwh),
        )
    for _, row in edges.iterrows():
        u = str(row.terminal_a)
        v = str(row.terminal_b)
        if u == v:
            continue
        length_m = float(row.length_m)
        if graph.has_edge(u, v) and graph.edges[u, v]["length_m"] <= length_m:
            continue
        graph.add_edge(
            u,
            v,
            edge_class=str(row.edge_class),
            length_m=length_m,
            length=length_m,
            edge_weight=float(row.edge_demand_kw),
            edge_demand_kw=float(row.edge_demand_kw),
            edge_transformer_count=int(row.edge_transformer_count),
            edge_size_kva=float(row.edge_size_kva),
            edge_consumer_count=float(row.edge_consumer_count),
            edge_yearly_kwh=float(row.edge_yearly_kwh),
            is_tie=False,
        )
    return graph


def read_final_network_tables(gpkg_path: Path) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]:
    backbone_edges = gpd.read_file(gpkg_path, layer="final_backbone_edges")
    tree_edges = gpd.read_file(gpkg_path, layer="final_tree_attachment_edges")
    transformer_nodes = gpd.read_file(gpkg_path, layer="final_transformer_nodes")
    source_nodes = gpd.read_file(gpkg_path, layer="final_source_nodes")
    edges = gpd.GeoDataFrame(
        list(backbone_edges.to_dict("records")) + list(tree_edges.to_dict("records")),
        geometry="geometry",
        crs=backbone_edges.crs,
    )
    nodes = transformer_nodes.copy()
    source_only = source_nodes[~source_nodes["terminal_id"].isin(nodes["terminal_id"])].copy()
    if not source_only.empty:
        nodes = gpd.GeoDataFrame(
            list(nodes.to_dict("records")) + list(source_only.to_dict("records")),
            geometry="geometry",
            crs=transformer_nodes.crs,
        )
    try:
        branch_nodes = gpd.read_file(gpkg_path, layer="final_street_branch_nodes")
    except Exception:
        branch_nodes = gpd.GeoDataFrame(geometry="geometry", crs=transformer_nodes.crs)
    if not branch_nodes.empty:
        nodes = gpd.GeoDataFrame(
            list(nodes.to_dict("records")) + list(branch_nodes.to_dict("records")),
            geometry="geometry",
            crs=transformer_nodes.crs,
        )
    return nodes, edges, backbone_edges, tree_edges


def _tree_terminal_row(
    *,
    road_node: str,
    transformer_by_road: dict[str, dict],
    node_pos: dict[str, tuple[float, float]],
    nearest_backbone_terminal: str | None = None,
    nearest_backbone_road_node: str | None = None,
    distance_to_backbone_m: float | None = None,
) -> dict:
    data = transformer_by_road[road_node]
    row = {
        "terminal_id": f"T:{road_node}",
        "road_node": road_node,
        "kind": "transformer",
        "terminal_class": "tree_terminal",
        "transformer_count": int(data["transformer_count"]),
        "source_count": 0,
        "size_kva": float(data["size_kva"]),
        "demand_kw": float(data["demand_kw"]),
        "consumer_count": float(data["consumer_count"]),
        "yearly_kwh": float(data["yearly_kwh"]),
        "geometry": line_from_road_nodes([road_node], node_pos).centroid,
    }
    if nearest_backbone_terminal is not None:
        row["nearest_backbone_terminal"] = nearest_backbone_terminal
    if nearest_backbone_road_node is not None:
        row["nearest_backbone_road_node"] = nearest_backbone_road_node
    if distance_to_backbone_m is not None:
        row["distance_to_backbone_m"] = float(distance_to_backbone_m)
    return row


def _star_tree_rows(
    *,
    nonbackbone_road_nodes: list[str],
    transformer_by_road: dict[str, dict],
    terminal_by_road: dict[str, str],
    dist: dict[str, float],
    paths: dict[str, list[str]],
    node_pos: dict[str, tuple[float, float]],
) -> tuple[list[dict], list[dict], list[dict], set[tuple[str, str]], int]:
    tree_node_rows = []
    tree_edge_rows = []
    tree_road_edges: set[tuple[str, str]] = set()
    unattached_tree_terminal_count = 0
    for road_node in nonbackbone_road_nodes:
        if road_node not in paths:
            unattached_tree_terminal_count += 1
            continue
        owner_road_node = clean(paths[road_node][0])
        owner_terminal = terminal_by_road.get(owner_road_node)
        if owner_terminal is None:
            unattached_tree_terminal_count += 1
            continue
        path = [clean(node) for node in paths[road_node]]
        for u, v in zip(path[:-1], path[1:]):
            tree_road_edges.add(edge_key(u, v))
        terminal_id = f"T:{road_node}"
        tree_node_rows.append(
            _tree_terminal_row(
                road_node=road_node,
                transformer_by_road=transformer_by_road,
                node_pos=node_pos,
                nearest_backbone_terminal=owner_terminal,
                nearest_backbone_road_node=owner_road_node,
                distance_to_backbone_m=float(dist[road_node]),
            )
        )
        tree_edge_rows.append(
            _tree_edge_row(
                terminal_a=terminal_id,
                terminal_b=owner_terminal,
                road_node_sequence=[road_node, owner_road_node],
                terminal_kind_a="transformer",
                terminal_kind_b="source_or_backbone",
                length_m=float(dist[road_node]),
                node_pos=node_pos,
            )
        )
    return tree_node_rows, [], tree_edge_rows, tree_road_edges, unattached_tree_terminal_count


def _street_forest_tree_rows(
    *,
    nonbackbone_road_nodes: list[str],
    transformer_by_road: dict[str, dict],
    terminal_by_road: dict[str, str],
    retained_backbone_road_nodes: set[str],
    road_graph: nx.Graph,
    dist: dict[str, float],
    paths: dict[str, list[str]],
    node_pos: dict[str, tuple[float, float]],
) -> tuple[list[dict], list[dict], list[dict], set[tuple[str, str]], int]:
    attached_transformer_nodes = []
    tree_road_edges: set[tuple[str, str]] = set()
    unattached_tree_terminal_count = 0
    for road_node in nonbackbone_road_nodes:
        if road_node not in paths:
            unattached_tree_terminal_count += 1
            continue
        attached_transformer_nodes.append(road_node)
        path = [clean(node) for node in paths[road_node]]
        for u, v in zip(path[:-1], path[1:]):
            tree_road_edges.add(edge_key(u, v))

    forest = nx.Graph()
    for u, v in tree_road_edges:
        if road_graph.has_edge(u, v):
            forest.add_edge(u, v, length_m=float(road_graph.edges[u, v]["length_m"]))

    tree_node_rows = []
    tree_terminal_by_road: dict[str, str] = {}
    for road_node in attached_transformer_nodes:
        path = [clean(node) for node in paths[road_node]]
        owner_road_node = clean(path[0])
        owner_terminal = terminal_by_road.get(owner_road_node)
        tree_terminal_by_road[road_node] = f"T:{road_node}"
        tree_node_rows.append(
            _tree_terminal_row(
                road_node=road_node,
                transformer_by_road=transformer_by_road,
                node_pos=node_pos,
                nearest_backbone_terminal=owner_terminal,
                nearest_backbone_road_node=owner_road_node,
                distance_to_backbone_m=float(dist[road_node]),
            )
        )

    boundary_roads = set(retained_backbone_road_nodes) | set(attached_transformer_nodes)
    boundary_roads.update(node for node, degree in forest.degree if degree != 2)

    branch_rows = []
    branch_terminal_by_road: dict[str, str] = {}
    for road_node in sorted(boundary_roads, key=str):
        if road_node in terminal_by_road or road_node in tree_terminal_by_road:
            continue
        if road_node not in forest:
            continue
        terminal_id = f"J:{road_node}"
        branch_terminal_by_road[road_node] = terminal_id
        branch_rows.append(
            {
                "terminal_id": terminal_id,
                "road_node": road_node,
                "kind": "street_branch",
                "terminal_class": "street_branch",
                "transformer_count": 0,
                "source_count": 0,
                "size_kva": 0.0,
                "demand_kw": 0.0,
                "consumer_count": 0.0,
                "yearly_kwh": 0.0,
                "geometry": line_from_road_nodes([road_node], node_pos).centroid,
            }
        )

    boundary_terminal_by_road = dict(terminal_by_road)
    boundary_terminal_by_road.update(tree_terminal_by_road)
    boundary_terminal_by_road.update(branch_terminal_by_road)
    tree_edge_rows = _contract_forest_edges(
        forest=forest,
        boundary_roads=boundary_roads,
        boundary_terminal_by_road=boundary_terminal_by_road,
        node_pos=node_pos,
    )
    return tree_node_rows, branch_rows, tree_edge_rows, tree_road_edges, unattached_tree_terminal_count


def _contract_forest_edges(
    *,
    forest: nx.Graph,
    boundary_roads: set[str],
    boundary_terminal_by_road: dict[str, str],
    node_pos: dict[str, tuple[float, float]],
) -> list[dict]:
    visited_edges: set[tuple[str, str]] = set()
    rows = []

    for start in sorted(boundary_roads & set(forest.nodes), key=str):
        for neighbor in sorted(forest.neighbors(start), key=str):
            ekey = edge_key(start, neighbor)
            if ekey in visited_edges:
                continue
            path, length_m = _walk_forest_chain(forest, start, neighbor, boundary_roads, visited_edges)
            end = path[-1]
            terminal_a = boundary_terminal_by_road.get(start)
            terminal_b = boundary_terminal_by_road.get(end)
            if terminal_a is None or terminal_b is None or terminal_a == terminal_b:
                continue
            rows.append(
                _tree_edge_row(
                    terminal_a=terminal_a,
                    terminal_b=terminal_b,
                    road_node_sequence=path,
                    terminal_kind_a=_terminal_kind_from_id(terminal_a),
                    terminal_kind_b=_terminal_kind_from_id(terminal_b),
                    length_m=length_m,
                    node_pos=node_pos,
                )
            )
    return rows


def _walk_forest_chain(
    forest: nx.Graph,
    start: str,
    neighbor: str,
    boundary_roads: set[str],
    visited_edges: set[tuple[str, str]],
) -> tuple[list[str], float]:
    path = [start, neighbor]
    visited_edges.add(edge_key(start, neighbor))
    total_length = float(forest.edges[start, neighbor]["length_m"])
    prev = start
    current = neighbor
    while current not in boundary_roads:
        candidates = [
            nxt
            for nxt in forest.neighbors(current)
            if nxt != prev and edge_key(current, nxt) not in visited_edges
        ]
        if not candidates:
            break
        nxt = sorted(candidates, key=str)[0]
        visited_edges.add(edge_key(current, nxt))
        total_length += float(forest.edges[current, nxt]["length_m"])
        path.append(nxt)
        prev, current = current, nxt
    return path, total_length


def _terminal_kind_from_id(terminal_id: str) -> str:
    if terminal_id.startswith("S:"):
        return "source"
    if terminal_id.startswith("J:"):
        return "street_branch"
    return "transformer"


def _tree_edge_row(
    *,
    terminal_a: str,
    terminal_b: str,
    road_node_sequence: list[str],
    terminal_kind_a: str,
    terminal_kind_b: str,
    length_m: float,
    node_pos: dict[str, tuple[float, float]],
) -> dict:
    return {
        "terminal_a": terminal_a,
        "terminal_b": terminal_b,
        "road_node_a": road_node_sequence[0],
        "road_node_b": road_node_sequence[-1],
        "terminal_kind_a": terminal_kind_a,
        "terminal_kind_b": terminal_kind_b,
        "edge_class": "tree_attachment",
        "length_m": float(length_m),
        "edge_transformer_count": 0,
        "edge_size_kva": 0.0,
        "edge_demand_kw": 0.0,
        "edge_consumer_count": 0.0,
        "edge_yearly_kwh": 0.0,
        "road_node_count": int(len(road_node_sequence)),
        "road_node_sequence": "|".join(road_node_sequence),
        "geometry": line_from_road_nodes(road_node_sequence, node_pos),
    }


def build_final_network_tables(
    backbone_summary: Path = BACKBONE_SUMMARY,
    *,
    tree_mode: str = "street_forest",
) -> dict:
    if tree_mode not in {"street_forest", "star"}:
        raise ValueError("tree_mode must be 'street_forest' or 'star'")

    backbone_gpkg, backbone_edges, backbone_transformer_nodes, backbone_source_nodes = load_backbone_solution(backbone_summary)
    branches = read_p2u_layer("StreetMap_branches")
    street_nodes = read_p2u_layer("StreetMap_nodes").copy()
    street_nodes["Node"] = street_nodes["Node"].map(clean)
    road_graph, node_pos = build_road_graph(branches, street_nodes)
    transformer_by_road, transformers, lv_assignment = transformers_by_street_node(street_nodes)

    retained_backbone_road_nodes = set(backbone_transformer_nodes["road_node"].astype(str).map(clean))
    retained_backbone_road_nodes.update(backbone_source_nodes["road_node"].astype(str).map(clean))
    retained_backbone_road_nodes.discard(SUPER_SOURCE)
    all_backbone_road_nodes = set(retained_backbone_road_nodes)
    for sequence in backbone_edges["terminal_sequence"].astype(str):
        all_backbone_road_nodes.update(sequence_road_nodes(sequence))

    missing_roads = sorted(node for node in retained_backbone_road_nodes if node not in road_graph)
    if missing_roads:
        raise ValueError(f"{len(missing_roads)} retained backbone road nodes are missing from the street graph")

    source_rows = solution_node_rows(backbone_source_nodes, "source", transformer_by_road, node_pos)
    backbone_node_rows = solution_node_rows(backbone_transformer_nodes, "backbone_boundary", transformer_by_road, node_pos)

    terminal_by_road = {
        clean(row.road_node): str(row.terminal_id)
        for _, row in backbone_transformer_nodes.iterrows()
    }
    terminal_by_road.update(
        {
            clean(row.road_node): str(row.terminal_id)
            for _, row in backbone_source_nodes.iterrows()
        }
    )

    dist, paths = nx.multi_source_dijkstra(road_graph, sources=list(retained_backbone_road_nodes), weight="length_m")

    nonbackbone_road_nodes = sorted(set(transformer_by_road) - all_backbone_road_nodes)
    if tree_mode == "star":
        tree_node_rows, tree_branch_rows, tree_edge_rows, tree_road_edges, unattached_tree_terminal_count = _star_tree_rows(
            nonbackbone_road_nodes=nonbackbone_road_nodes,
            transformer_by_road=transformer_by_road,
            terminal_by_road=terminal_by_road,
            dist=dist,
            paths=paths,
            node_pos=node_pos,
        )
    else:
        tree_node_rows, tree_branch_rows, tree_edge_rows, tree_road_edges, unattached_tree_terminal_count = _street_forest_tree_rows(
            nonbackbone_road_nodes=nonbackbone_road_nodes,
            transformer_by_road=transformer_by_road,
            terminal_by_road=terminal_by_road,
            retained_backbone_road_nodes=retained_backbone_road_nodes,
            road_graph=road_graph,
            dist=dist,
            paths=paths,
            node_pos=node_pos,
        )

    final_node_rows = backbone_node_rows + source_rows + tree_node_rows + tree_branch_rows
    final_nodes = gpd.GeoDataFrame(final_node_rows, geometry="geometry", crs=SOURCE_CRS)
    final_transformers = final_nodes[final_nodes["kind"].isin(["transformer", "source_transformer"])].copy()
    final_sources = final_nodes[final_nodes["source_count"].astype(float) > 0].copy()
    final_branch_nodes = final_nodes[final_nodes["kind"] == "street_branch"].copy()

    backbone_edge_rows = []
    for _, row in backbone_edges.iterrows():
        loads = edge_internal_load(str(row.terminal_sequence), transformer_by_road)
        backbone_edge_rows.append(
            {
                "terminal_a": str(row.terminal_a),
                "terminal_b": str(row.terminal_b),
                "road_node_a": clean(row.road_node_a),
                "road_node_b": clean(row.road_node_b),
                "terminal_kind_a": str(row.terminal_kind_a),
                "terminal_kind_b": str(row.terminal_kind_b),
                "edge_class": "backbone",
                "length_m": float(row.length_m),
                "terminal_sequence": str(row.terminal_sequence),
                **loads,
                "geometry": line_from_road_nodes([clean(row.road_node_a), clean(row.road_node_b)], node_pos),
            }
        )
    final_backbone_edges = gpd.GeoDataFrame(backbone_edge_rows, geometry="geometry", crs=SOURCE_CRS)
    final_tree_edges = gpd.GeoDataFrame(tree_edge_rows, geometry="geometry", crs=SOURCE_CRS)
    final_edges = gpd.GeoDataFrame(backbone_edge_rows + tree_edge_rows, geometry="geometry", crs=SOURCE_CRS)

    tree_road_union_length = float(
        sum(road_graph.edges[u, v]["length_m"] for u, v in tree_road_edges if road_graph.has_edge(u, v))
    )
    metadata = {
        "input_backbone_gpkg": str(backbone_gpkg),
        "final_graph_definition": (
            "switch-section graph: graph nodes are transformer/source terminals "
            "plus optional zero-load street branch nodes; "
            "degree-2 transformer chains are stored as edge load; "
            f"tree mode is {tree_mode}"
        ),
        "tree_mode": tree_mode,
        "retained_backbone_road_nodes": int(len(retained_backbone_road_nodes)),
        "all_backbone_terminal_sequence_road_nodes": int(len(all_backbone_road_nodes)),
        "total_transformers_in_data": int(len(transformers)),
        "unattached_tree_terminal_count": int(unattached_tree_terminal_count),
        "tree_attachment_physical_road_union_length_m": tree_road_union_length,
        "street_branch_nodes": int(len(final_branch_nodes)),
        "lv_assignment": lv_assignment,
    }
    return {
        "final_nodes": final_nodes,
        "final_transformer_nodes": final_transformers,
        "final_source_nodes": final_sources,
        "final_street_branch_nodes": final_branch_nodes,
        "final_edges": final_edges,
        "final_backbone_edges": final_backbone_edges,
        "final_tree_attachment_edges": final_tree_edges,
        "metadata": metadata,
    }


def write_final_network(
    tables: dict,
    *,
    output_gpkg: Path = FINAL_NETWORK_GPKG,
    metadata_json: Path = FINAL_NETWORK_METADATA,
) -> dict:
    output_gpkg.parent.mkdir(parents=True, exist_ok=True)
    temp_gpkg = output_gpkg.with_name(f"{output_gpkg.stem}.tmp{output_gpkg.suffix}")
    if temp_gpkg.exists():
        temp_gpkg.unlink()

    layers = {
        "final_backbone_edges": tables["final_backbone_edges"],
        "final_tree_attachment_edges": tables["final_tree_attachment_edges"],
        "final_transformer_nodes": tables["final_transformer_nodes"],
        "final_source_nodes": tables["final_source_nodes"],
        "final_street_branch_nodes": tables.get("final_street_branch_nodes", gpd.GeoDataFrame()),
    }
    for layer, gdf in layers.items():
        if not gdf.empty:
            gdf.to_crs(TARGET_CRS).to_file(temp_gpkg, layer=layer, driver="GPKG")

    actual_gpkg = output_gpkg
    try:
        if output_gpkg.exists():
            output_gpkg.unlink()
        temp_gpkg.replace(output_gpkg)
    except PermissionError:
        actual_gpkg = output_gpkg.with_name(f"{output_gpkg.stem}_{int(time.time())}{output_gpkg.suffix}")
        temp_gpkg.replace(actual_gpkg)

    metadata = dict(tables["metadata"])
    metadata["output_gpkg"] = str(actual_gpkg)
    metadata["metadata_json"] = str(metadata_json)
    metadata_json.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def build_and_write_final_network(
    *,
    backbone_summary: Path = BACKBONE_SUMMARY,
    output_gpkg: Path = FINAL_NETWORK_GPKG,
    metadata_json: Path = FINAL_NETWORK_METADATA,
    tree_mode: str = "street_forest",
) -> dict:
    tables = build_final_network_tables(backbone_summary, tree_mode=tree_mode)
    return write_final_network(tables, output_gpkg=output_gpkg, metadata_json=metadata_json)
