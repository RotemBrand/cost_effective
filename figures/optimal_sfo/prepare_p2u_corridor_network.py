from __future__ import annotations

import heapq
import html
import json
from collections import Counter, defaultdict
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import LineString


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "power" / "better_grids" / "SFO" / "P2U"
OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo"
GPKG_PATH = OUTPUT_DIR / "p2u_terminal_corridors_road2_k10_3857.gpkg"
QGS_PATH = OUTPUT_DIR / "p2u_terminal_corridors_road2_k10.qgs"
SUMMARY_PATH = OUTPUT_DIR / "p2u_terminal_corridors_road2_k10_summary.md"
SUMMARY_JSON_PATH = OUTPUT_DIR / "p2u_terminal_corridors_road2_k10_summary.json"

SOURCE_CRS = "EPSG:32610"
TARGET_CRS = "EPSG:3857"
SUPER_SOURCE = "__SOURCE__"
K_NEAREST_TERMINAL_EDGES = 10


def _clean(value) -> str:
    return str(value).strip()


def _read_layer(name: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(DATA_DIR / f"{name}.shp")
    if gdf.crs is None:
        gdf = gdf.set_crs(SOURCE_CRS)
    return gdf


def _nearest_nodes(points: gpd.GeoDataFrame, street_nodes: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    points = points.to_crs(street_nodes.crs)
    point_geom = points.geometry.centroid
    node_geom = street_nodes.geometry
    tree = cKDTree(np.column_stack([node_geom.x.to_numpy(), node_geom.y.to_numpy()]))
    distances, indices = tree.query(np.column_stack([point_geom.x.to_numpy(), point_geom.y.to_numpy()]), k=1)
    return street_nodes.iloc[indices]["Node"].to_numpy(), distances


def _build_road_graph(
    branches: gpd.GeoDataFrame, street_nodes: gpd.GeoDataFrame
) -> tuple[nx.Graph, dict[str, tuple[float, float]]]:
    node_pos = {
        _clean(row.Node): (float(row.geometry.x), float(row.geometry.y))
        for _, row in street_nodes.iterrows()
    }
    graph = nx.Graph()
    duplicates = 0
    for _, row in branches.iterrows():
        u = _clean(row.Node_A)
        v = _clean(row.Node_B)
        if u == v:
            continue
        length_m = float(row.geometry.length)
        if graph.has_edge(u, v):
            duplicates += 1
            if length_m < graph.edges[u, v]["length_m"]:
                graph.edges[u, v]["length_m"] = length_m
        else:
            graph.add_edge(u, v, length_m=length_m)
    graph.graph["duplicate_edges"] = duplicates
    return graph, node_pos


def _contract_sources(graph: nx.Graph, source_nodes: set[str]) -> nx.Graph:
    contracted = nx.Graph()
    for u, v, data in graph.edges(data=True):
        cu = SUPER_SOURCE if u in source_nodes else u
        cv = SUPER_SOURCE if v in source_nodes else v
        if cu == cv:
            continue
        length_m = float(data.get("length_m", 1.0))
        if contracted.has_edge(cu, cv):
            contracted.edges[cu, cv]["length_m"] = min(
                contracted.edges[cu, cv]["length_m"], length_m
            )
        else:
            contracted.add_edge(cu, cv, length_m=length_m)
    return contracted


def _source_side_2edge_original_nodes(
    road_graph: nx.Graph, source_nodes: set[str]
) -> set[str]:
    contracted = _contract_sources(road_graph, source_nodes)
    source_blocks = [
        set(block)
        for block in nx.k_edge_components(contracted, k=2)
        if SUPER_SOURCE in set(block)
    ]
    if not source_blocks:
        raise RuntimeError("No source-side 2-edge-connected street component was found.")
    source_block = max(source_blocks, key=len)
    selected = set(source_block)
    selected.discard(SUPER_SOURCE)
    if SUPER_SOURCE in source_block:
        selected.update(source_nodes)
    return {node for node in selected if node in road_graph}


def _terminal_id(kind: str, road_node: str) -> str:
    return f"{kind}:{road_node}"


def _build_terminal_tables(
    street_nodes: gpd.GeoDataFrame,
    transformers: gpd.GeoDataFrame,
    substations: gpd.GeoDataFrame,
    selected_nodes: set[str],
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, dict[str, str], dict[str, dict]]:
    transformer_street_nodes, transformer_distances = _nearest_nodes(transformers, street_nodes)
    source_street_nodes, source_distances = _nearest_nodes(substations, street_nodes)

    source_nodes = set(map(_clean, source_street_nodes))
    transformer_rows = transformers.copy()
    transformer_rows["transformer"] = transformer_rows["Node"].map(_clean)
    transformer_rows["street_node"] = list(map(_clean, transformer_street_nodes))
    transformer_rows["nearest_street_distance_m"] = transformer_distances.astype(float)
    transformer_rows["size_kva"] = transformer_rows["Size_kVA"].astype(float)
    transformer_rows["in_selected_2edge_component"] = transformer_rows["street_node"].isin(selected_nodes)

    source_rows = substations.copy()
    source_rows["source"] = source_rows["Node"].map(_clean)
    source_rows["street_node"] = list(map(_clean, source_street_nodes))
    source_rows["nearest_street_distance_m"] = source_distances.astype(float)
    source_rows["in_selected_2edge_component"] = source_rows["street_node"].isin(selected_nodes)

    terminal_by_node_data: dict[str, dict] = defaultdict(
        lambda: {
            "transformer_count": 0,
            "source_count": 0,
            "size_kva": 0.0,
        }
    )

    for road_node, group in transformer_rows[transformer_rows["in_selected_2edge_component"]].groupby("street_node"):
        terminal_by_node_data[road_node]["transformer_count"] += int(len(group))
        terminal_by_node_data[road_node]["size_kva"] += float(group["size_kva"].sum())

    for road_node, group in source_rows[source_rows["in_selected_2edge_component"]].groupby("street_node"):
        terminal_by_node_data[road_node]["source_count"] += int(len(group))

    terminal_by_road_node: dict[str, str] = {}
    terminal_data: dict[str, dict] = {}
    for road_node, data in terminal_by_node_data.items():
        has_transformer = data["transformer_count"] > 0
        has_source = data["source_count"] > 0
        if has_transformer and has_source:
            kind = "source_transformer"
            prefix = "ST"
        elif has_source:
            kind = "source"
            prefix = "S"
        else:
            kind = "transformer"
            prefix = "T"
        terminal = _terminal_id(prefix, road_node)
        terminal_by_road_node[road_node] = terminal
        terminal_data[terminal] = {
            "terminal_id": terminal,
            "kind": kind,
            "road_node": road_node,
            "transformer_count": int(data["transformer_count"]),
            "source_count": int(data["source_count"]),
            "size_kva": float(data["size_kva"]),
        }

    return transformer_rows, source_rows, terminal_by_road_node, terminal_data


def _multi_source_terminal_voronoi(
    graph: nx.Graph, terminal_by_road_node: dict[str, str]
) -> tuple[dict[str, str], dict[str, float], dict[str, str | None]]:
    owner: dict[str, str] = {}
    dist: dict[str, float] = {}
    parent: dict[str, str | None] = {}
    queue: list[tuple[float, str, str]] = []

    for road_node, terminal in terminal_by_road_node.items():
        owner[road_node] = terminal
        dist[road_node] = 0.0
        parent[road_node] = None
        heapq.heappush(queue, (0.0, terminal, road_node))

    while queue:
        distance, terminal, node = heapq.heappop(queue)
        if distance != dist.get(node) or terminal != owner.get(node):
            continue
        for neighbor, edge_data in graph[node].items():
            candidate = distance + float(edge_data.get("length_m", 1.0))
            current = dist.get(neighbor)
            if current is None or candidate < current:
                dist[neighbor] = candidate
                owner[neighbor] = terminal
                parent[neighbor] = node
                heapq.heappush(queue, (candidate, terminal, neighbor))

    return owner, dist, parent


def _path_owner_to_node(node: str, parent: dict[str, str | None]) -> list[str]:
    path = [node]
    while parent[path[-1]] is not None:
        path.append(parent[path[-1]])
    path.reverse()
    return path


def _line_from_node_sequence(nodes: list[str], node_pos: dict[str, tuple[float, float]]) -> LineString:
    coords = []
    for node in nodes:
        if node not in node_pos:
            continue
        coord = node_pos[node]
        if not coords or coords[-1] != coord:
            coords.append(coord)
    if len(coords) == 1:
        coords.append(coords[0])
    return LineString(coords)


def _straight_terminal_line(
    terminal_a: str,
    terminal_b: str,
    terminal_data: dict[str, dict],
    node_pos: dict[str, tuple[float, float]],
) -> LineString:
    return _line_from_node_sequence(
        [terminal_data[terminal_a]["road_node"], terminal_data[terminal_b]["road_node"]],
        node_pos,
    )


def _add_or_update_terminal_edge(
    graph: nx.Graph,
    terminal_a: str,
    terminal_b: str,
    length_m: float,
    component_id: int,
    internal_road_node_count: int,
) -> None:
    if terminal_a == terminal_b:
        return
    if graph.has_edge(terminal_a, terminal_b) and graph.edges[terminal_a, terminal_b]["length_m"] <= length_m:
        return
    graph.add_edge(
        terminal_a,
        terminal_b,
        length_m=float(length_m),
        nonterminal_component_id=int(component_id),
        internal_road_node_count=int(internal_road_node_count),
    )


def _build_terminal_minor_graph(
    road_graph: nx.Graph,
    terminal_by_road_node: dict[str, str],
    terminal_data: dict[str, dict],
) -> tuple[nx.Graph, dict[str, int | float]]:
    """Build the graph whose nodes are terminals and whose edges skip non-terminals.

    An edge exists when its two endpoint terminals are connected by a road path
    whose internal nodes are not transformer/source terminals.
    """
    terminal_road_nodes = set(terminal_by_road_node)
    graph = nx.Graph()
    for terminal, data in terminal_data.items():
        graph.add_node(terminal, **data)

    direct_terminal_edges = 0
    for u, v, data in road_graph.edges(data=True):
        if u in terminal_road_nodes and v in terminal_road_nodes:
            direct_terminal_edges += 1
            _add_or_update_terminal_edge(
                graph,
                terminal_by_road_node[u],
                terminal_by_road_node[v],
                float(data.get("length_m", 1.0)),
                component_id=-1,
                internal_road_node_count=0,
            )

    nonterminal_nodes = set(road_graph.nodes) - terminal_road_nodes
    nonterminal_graph = road_graph.subgraph(nonterminal_nodes)
    component_count = 0
    boundary_histogram = Counter()
    max_boundary_terminals = 0

    for component_id, component_nodes_iter in enumerate(nx.connected_components(nonterminal_graph)):
        component_nodes = set(component_nodes_iter)
        boundary_adjacency: dict[str, list[tuple[str, float]]] = defaultdict(list)
        for node in component_nodes:
            for neighbor, edge_data in road_graph[node].items():
                if neighbor in terminal_by_road_node:
                    boundary_adjacency[terminal_by_road_node[neighbor]].append(
                        (node, float(edge_data.get("length_m", 1.0)))
                    )

        boundary_terminals = sorted(boundary_adjacency)
        if len(boundary_terminals) < 2:
            continue
        component_count += 1
        boundary_histogram[len(boundary_terminals)] += 1
        max_boundary_terminals = max(max_boundary_terminals, len(boundary_terminals))

        augmented = road_graph.subgraph(component_nodes).copy()
        for terminal, attachments in boundary_adjacency.items():
            augmented.add_node(terminal)
            for road_node, length_m in attachments:
                if augmented.has_edge(terminal, road_node):
                    augmented.edges[terminal, road_node]["length_m"] = min(
                        augmented.edges[terminal, road_node]["length_m"], length_m
                    )
                else:
                    augmented.add_edge(terminal, road_node, length_m=length_m)

        for i, terminal_a in enumerate(boundary_terminals):
            distances = nx.single_source_dijkstra_path_length(
                augmented, terminal_a, weight="length_m"
            )
            for terminal_b in boundary_terminals[i + 1 :]:
                if terminal_b not in distances:
                    continue
                _add_or_update_terminal_edge(
                    graph,
                    terminal_a,
                    terminal_b,
                    distances[terminal_b],
                    component_id=component_id,
                    internal_road_node_count=len(component_nodes),
                )

    stats = {
        "terminal_nodes_before_selection": graph.number_of_nodes(),
        "terminal_edges_before_selection": graph.number_of_edges(),
        "direct_terminal_road_edges": direct_terminal_edges,
        "nonterminal_components_with_2plus_boundary_terminals": component_count,
        "max_boundary_terminals_on_nonterminal_component": max_boundary_terminals,
        "boundary_terminal_histogram": dict(sorted(boundary_histogram.items())),
    }
    return graph, stats


def _source_side_terminal_graph(
    terminal_graph: nx.Graph,
) -> nx.Graph:
    source_like = {"source", "source_transformer"}
    source_terminals = {
        node for node, data in terminal_graph.nodes(data=True) if data["kind"] in source_like
    }
    contracted = nx.Graph()
    for u, v, data in terminal_graph.edges(data=True):
        cu = SUPER_SOURCE if u in source_terminals else u
        cv = SUPER_SOURCE if v in source_terminals else v
        if cu == cv:
            continue
        length_m = float(data.get("length_m", 1.0))
        if contracted.has_edge(cu, cv):
            contracted.edges[cu, cv]["length_m"] = min(contracted.edges[cu, cv]["length_m"], length_m)
        else:
            contracted.add_edge(cu, cv, length_m=length_m)

    source_blocks = [
        set(block)
        for block in nx.k_edge_components(contracted, k=2)
        if SUPER_SOURCE in set(block)
    ]
    if not source_blocks:
        raise RuntimeError("No source-side 2-edge-connected terminal component was found.")
    source_block = max(source_blocks, key=len)

    selected = nx.Graph()
    for u, v, data in terminal_graph.edges(data=True):
        cu = SUPER_SOURCE if u in source_terminals else u
        cv = SUPER_SOURCE if v in source_terminals else v
        if cu in source_block and cv in source_block:
            selected.add_node(u, **terminal_graph.nodes[u])
            selected.add_node(v, **terminal_graph.nodes[v])
            selected.add_edge(u, v, **data)
    return selected


def _prune_terminal_graph_k_nearest(graph: nx.Graph, k: int) -> nx.Graph:
    """Keep the symmetric k shortest candidate edges incident to each terminal."""
    keep_edges: set[tuple[str, str]] = set()

    def edge_key(u: str, v: str) -> tuple[str, str]:
        return tuple(sorted((u, v)))

    for node in graph.nodes:
        incident = sorted(
            graph.edges(node, data=True),
            key=lambda item: (float(item[2].get("length_m", 1.0)), str(item[1])),
        )
        for u, v, _ in incident[:k]:
            keep_edges.add(edge_key(u, v))

    pruned = nx.Graph()
    for u, data in graph.nodes(data=True):
        pruned.add_node(u, **data)
    for u, v in keep_edges:
        pruned.add_edge(u, v, **graph.edges[u, v])

    isolated = [node for node in pruned.nodes if pruned.degree(node) == 0]
    pruned.remove_nodes_from(isolated)
    return pruned


def _source_side_2edge_status(graph: nx.Graph) -> dict[str, int | bool]:
    source_like = {"source", "source_transformer"}
    source_nodes = {node for node, data in graph.nodes(data=True) if data["kind"] in source_like}
    selected = _source_side_terminal_graph(graph)
    return {
        "is_entire_pruned_graph_source_side_2edge": (
            selected.number_of_nodes() == graph.number_of_nodes()
            and selected.number_of_edges() == graph.number_of_edges()
        ),
        "source_side_2edge_nodes_after_pruning": selected.number_of_nodes(),
        "source_side_2edge_edges_after_pruning": selected.number_of_edges(),
        "pruned_terminal_nodes": graph.number_of_nodes(),
        "pruned_terminal_edges": graph.number_of_edges(),
        "source_terminals_after_pruning": len(source_nodes),
    }


def _terminal_nodes_gdf(
    graph: nx.Graph,
    terminal_data: dict[str, dict],
    node_pos: dict[str, tuple[float, float]],
) -> gpd.GeoDataFrame:
    rows = []
    for node in graph.nodes:
        data = terminal_data[node]
        road_node = data["road_node"]
        rows.append(
            {
                "terminal_id": node,
                "road_node": road_node,
                "kind": data["kind"],
                "transformer_count": int(data["transformer_count"]),
                "source_count": int(data["source_count"]),
                "size_kva": float(data["size_kva"]),
                "degree": int(graph.degree(node)),
                "geometry": _line_from_node_sequence([road_node], node_pos).centroid,
            }
        )
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=SOURCE_CRS)


def _terminal_edges_gdf(
    graph: nx.Graph,
    terminal_data: dict[str, dict],
    node_pos: dict[str, tuple[float, float]],
) -> gpd.GeoDataFrame:
    rows = []
    for u, v, data in graph.edges(data=True):
        rows.append(
            {
                "terminal_a": u,
                "terminal_b": v,
                "road_node_a": terminal_data[u]["road_node"],
                "road_node_b": terminal_data[v]["road_node"],
                "terminal_kind_a": terminal_data[u]["kind"],
                "terminal_kind_b": terminal_data[v]["kind"],
                "length_m": float(data.get("length_m", 1.0)),
                "nonterminal_component_id": int(data.get("nonterminal_component_id", -1)),
                "internal_road_node_count": int(data.get("internal_road_node_count", 0)),
                "geometry": _straight_terminal_line(u, v, terminal_data, node_pos),
            }
        )
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=SOURCE_CRS)


def _extract_terminal_corridors(
    graph: nx.Graph,
    node_pos: dict[str, tuple[float, float]],
    terminal_by_road_node: dict[str, str],
    terminal_data: dict[str, dict],
) -> gpd.GeoDataFrame:
    owner, dist, parent = _multi_source_terminal_voronoi(graph, terminal_by_road_node)
    best_by_pair: dict[tuple[str, str], dict] = {}

    for u, v, edge_data in graph.edges(data=True):
        owner_u = owner.get(u)
        owner_v = owner.get(v)
        if owner_u is None or owner_v is None or owner_u == owner_v:
            continue
        pair = tuple(sorted((owner_u, owner_v)))
        length_m = dist[u] + float(edge_data.get("length_m", 1.0)) + dist[v]
        current = best_by_pair.get(pair)
        if current is not None and length_m >= current["length_m"]:
            continue

        path_u = _path_owner_to_node(u, parent)
        path_v = _path_owner_to_node(v, parent)
        node_sequence = path_u + list(reversed(path_v))
        best_by_pair[pair] = {
            "terminal_a": pair[0],
            "terminal_b": pair[1],
            "road_node_a": terminal_data[pair[0]]["road_node"],
            "road_node_b": terminal_data[pair[1]]["road_node"],
            "terminal_kind_a": terminal_data[pair[0]]["kind"],
            "terminal_kind_b": terminal_data[pair[1]]["kind"],
            "transformer_count_a": terminal_data[pair[0]]["transformer_count"],
            "transformer_count_b": terminal_data[pair[1]]["transformer_count"],
            "size_kva_a": terminal_data[pair[0]]["size_kva"],
            "size_kva_b": terminal_data[pair[1]]["size_kva"],
            "length_m": float(length_m),
            "road_node_count": len(node_sequence),
            "geometry": _line_from_node_sequence(
                [terminal_data[pair[0]]["road_node"], terminal_data[pair[1]]["road_node"]],
                node_pos,
            ),
        }

    return gpd.GeoDataFrame(list(best_by_pair.values()), geometry="geometry", crs=SOURCE_CRS)


def _build_contracted_corridors(
    corridor_rows: gpd.GeoDataFrame,
    terminal_data: dict[str, dict],
    node_pos: dict[str, tuple[float, float]],
) -> gpd.GeoDataFrame:
    graph = nx.Graph()
    for terminal, data in terminal_data.items():
        graph.add_node(terminal, **data)
    for _, row in corridor_rows.iterrows():
        graph.add_edge(
            row["terminal_a"],
            row["terminal_b"],
            length_m=float(row["length_m"]),
            road_node_count=int(row.get("road_node_count", row.get("internal_road_node_count", 0))),
        )

    source_like = {"source", "source_transformer"}
    boundary_nodes = {
        node
        for node in graph.nodes
        if graph.nodes[node]["kind"] in source_like or graph.degree(node) != 2
    }
    if not boundary_nodes and graph.number_of_nodes() > 0:
        boundary_nodes.add(max(graph.nodes, key=lambda n: graph.nodes[n].get("size_kva", 0.0)))

    visited_edges: set[tuple[str, str]] = set()
    rows = []

    def edge_key(u: str, v: str) -> tuple[str, str]:
        return tuple(sorted((u, v)))

    def trace_chain(start: str, next_node: str) -> tuple[list[str], float, int]:
        sequence = [start, next_node]
        total_length = float(graph.edges[start, next_node]["length_m"])
        total_road_nodes = int(graph.edges[start, next_node].get("road_node_count", 2))
        prev, current = start, next_node
        visited_edges.add(edge_key(prev, current))

        while current not in boundary_nodes and graph.degree(current) == 2:
            neighbors = list(graph.neighbors(current))
            nxt = neighbors[0] if neighbors[1] == prev else neighbors[1]
            ekey = edge_key(current, nxt)
            if ekey in visited_edges:
                break
            visited_edges.add(ekey)
            total_length += float(graph.edges[current, nxt]["length_m"])
            total_road_nodes += max(0, int(graph.edges[current, nxt].get("road_node_count", 2)) - 1)
            sequence.append(nxt)
            prev, current = current, nxt
        return sequence, total_length, total_road_nodes

    for start in sorted(boundary_nodes):
        for neighbor in graph.neighbors(start):
            if edge_key(start, neighbor) in visited_edges:
                continue
            sequence, total_length, total_road_nodes = trace_chain(start, neighbor)
            end = sequence[-1]
            internal = sequence[1:-1]
            road_a = graph.nodes[start]["road_node"]
            road_b = graph.nodes[end]["road_node"]
            rows.append(
                {
                    "terminal_a": start,
                    "terminal_b": end,
                    "road_node_a": road_a,
                    "road_node_b": road_b,
                    "terminal_kind_a": graph.nodes[start]["kind"],
                    "terminal_kind_b": graph.nodes[end]["kind"],
                    "edge_transformer_count": int(sum(graph.nodes[n].get("transformer_count", 0) for n in internal)),
                    "edge_source_count": int(sum(graph.nodes[n].get("source_count", 0) for n in internal)),
                    "edge_size_kva": float(sum(graph.nodes[n].get("size_kva", 0.0) for n in internal)),
                    "chain_terminal_count": len(sequence),
                    "internal_terminal_count": len(internal),
                    "length_m": float(total_length),
                    "road_node_count": int(total_road_nodes),
                    "terminal_sequence": "|".join(sequence),
                    "geometry": _line_from_node_sequence([road_a, road_b], node_pos),
                }
            )

    # Handle any remaining all-degree-2 components.
    for u, v in graph.edges:
        if edge_key(u, v) in visited_edges:
            continue
        sequence, total_length, total_road_nodes = trace_chain(u, v)
        end = sequence[-1]
        internal = sequence[1:-1]
        road_a = graph.nodes[u]["road_node"]
        road_b = graph.nodes[end]["road_node"]
        rows.append(
            {
                "terminal_a": u,
                "terminal_b": end,
                "road_node_a": road_a,
                "road_node_b": road_b,
                "terminal_kind_a": graph.nodes[u]["kind"],
                "terminal_kind_b": graph.nodes[end]["kind"],
                "edge_transformer_count": int(sum(graph.nodes[n].get("transformer_count", 0) for n in internal)),
                "edge_source_count": int(sum(graph.nodes[n].get("source_count", 0) for n in internal)),
                "edge_size_kva": float(sum(graph.nodes[n].get("size_kva", 0.0) for n in internal)),
                "chain_terminal_count": len(sequence),
                "internal_terminal_count": len(internal),
                "length_m": float(total_length),
                "road_node_count": int(total_road_nodes),
                "terminal_sequence": "|".join(sequence),
                "geometry": _line_from_node_sequence([road_a, road_b], node_pos),
            }
        )

    return gpd.GeoDataFrame(rows, geometry="geometry", crs=SOURCE_CRS)


def _write_layer(gdf: gpd.GeoDataFrame, layer: str) -> None:
    if gdf.empty:
        return
    gdf.to_crs(TARGET_CRS).to_file(GPKG_PATH, layer=layer, driver="GPKG")


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
                <Option name="outline_width" type="QString" value="0.35"/>
                <Option name="outline_width_unit" type="QString" value="MM"/>
                <Option name="size" type="QString" value="{size}"/>
                <Option name="size_unit" type="QString" value="MM"/>
              </Option>
            </layer>
          </symbol>
        </symbols>
      </renderer-v2>"""


def _map_layer(layer: str, name: str, layer_id: str, geometry: str, renderer: str) -> str:
    return f"""  <maplayer type="vector" geometry="{geometry}">
    <id>{layer_id}</id>
    <datasource>./{GPKG_PATH.name}|layername={layer}</datasource>
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


def _write_qgs() -> None:
    layers = [
        ("ilp_edges", "ILP edges", "ilp_edges_5f76b9cd", "Line", _line_symbol("163,76,156,255", "0.65"), "1"),
        ("ilp_transformer_nodes", "ILP transformer nodes", "ilp_transformer_nodes_8755a858", "Point", _point_symbol("255,255,255,0", "203,64,60,255", "1.25", False), "1"),
        ("ilp_source_nodes", "ILP source nodes", "ilp_source_nodes_1b7c0c88", "Point", _point_symbol("18,18,18,255", "255,255,255,255", "2.5", True), "1"),
        ("raw_terminal_edges", "Raw terminal-minor edges", "raw_terminal_edges_2988614c", "Line", _line_symbol("80,120,190,180", "0.35"), "0"),
        ("raw_terminal_nodes", "Raw terminal-minor nodes", "raw_terminal_nodes_d0e94731", "Point", _point_symbol("255,255,255,0", "120,120,120,180", "0.9", False), "0"),
    ]
    tree = "\n".join(
        f'      <layer-tree-layer checked="{checked}" id="{layer_id}" name="{html.escape(name)}" source="./{GPKG_PATH.name}|layername={layer}"/>'
        for layer, name, layer_id, _, _, checked in layers
    )
    project_layers = "\n".join(
        _map_layer(layer, name, layer_id, geometry, renderer)
        for layer, name, layer_id, geometry, renderer, _ in layers
    )
    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="P2U terminal corridors">
  <homePath path="."/>
  <title>P2U terminal corridors</title>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="P2U terminal corridor graph">
{tree}
    </layer-tree-group>
  </layer-tree-group>
  <projectlayers>
{project_layers}
  </projectlayers>
</qgis>
"""
    QGS_PATH.write_text(qgs, encoding="utf-8")


def build_corridor_outputs() -> dict[str, int | float]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if GPKG_PATH.exists():
        GPKG_PATH.unlink()

    branches = _read_layer("StreetMap_branches")
    street_nodes = _read_layer("StreetMap_nodes").copy()
    transformers = _read_layer("DistribTransf_N")
    substations = _read_layer("HVMVSubstation_N")
    street_nodes["Node"] = street_nodes["Node"].map(_clean)

    road_graph, node_pos = _build_road_graph(branches, street_nodes)
    source_street_nodes, _ = _nearest_nodes(substations, street_nodes)
    selected_road_nodes = _source_side_2edge_original_nodes(
        road_graph, set(map(_clean, source_street_nodes))
    )
    selected_road_graph = road_graph.subgraph(selected_road_nodes).copy()
    transformer_rows, source_rows, terminal_by_road_node, terminal_data = _build_terminal_tables(
        street_nodes, transformers, substations, selected_road_nodes
    )

    terminal_graph, terminal_minor_stats = _build_terminal_minor_graph(
        selected_road_graph, terminal_by_road_node, terminal_data
    )
    pruned_terminal_graph = _prune_terminal_graph_k_nearest(
        terminal_graph, K_NEAREST_TERMINAL_EDGES
    )
    pruning_status = _source_side_2edge_status(pruned_terminal_graph)
    selected_terminal_graph = _source_side_terminal_graph(pruned_terminal_graph)
    raw_terminal_edges = _terminal_edges_gdf(selected_terminal_graph, terminal_data, node_pos)
    raw_terminal_nodes = _terminal_nodes_gdf(selected_terminal_graph, terminal_data, node_pos)
    ilp_edges = _build_contracted_corridors(raw_terminal_edges, terminal_data, node_pos)

    ilp_node_ids = set(ilp_edges["terminal_a"]) | set(ilp_edges["terminal_b"]) if len(ilp_edges) else set()
    ilp_node_graph = selected_terminal_graph.subgraph(ilp_node_ids).copy()
    ilp_nodes = _terminal_nodes_gdf(ilp_node_graph, terminal_data, node_pos)
    ilp_transformer_nodes = ilp_nodes[ilp_nodes["kind"].isin(["transformer", "source_transformer"])].copy()
    ilp_source_nodes = ilp_nodes[ilp_nodes["kind"].isin(["source", "source_transformer"])].copy()

    _write_layer(ilp_edges, "ilp_edges")
    _write_layer(ilp_transformer_nodes, "ilp_transformer_nodes")
    _write_layer(ilp_source_nodes, "ilp_source_nodes")
    _write_layer(raw_terminal_edges, "raw_terminal_edges")
    _write_layer(raw_terminal_nodes, "raw_terminal_nodes")
    _write_qgs()

    kind_counts = Counter()
    for _, row in raw_terminal_edges.iterrows():
        kinds = tuple(sorted((row["terminal_kind_a"], row["terminal_kind_b"])))
        kind_counts["-".join(kinds)] += 1

    source_like = {"source", "source_transformer"}
    transformer_like = {"transformer", "source_transformer"}
    source_touching_edges = int(
        sum(
            1
            for _, row in raw_terminal_edges.iterrows()
            if row["terminal_kind_a"] in source_like or row["terminal_kind_b"] in source_like
        )
    )
    transformer_only_edges = int(
        sum(
            1
            for _, row in raw_terminal_edges.iterrows()
            if row["terminal_kind_a"] == "transformer" and row["terminal_kind_b"] == "transformer"
        )
    )

    summary = {
        "road_graph_nodes": road_graph.number_of_nodes(),
        "road_graph_edges": road_graph.number_of_edges(),
        "selected_source_side_road_2edge_nodes": selected_road_graph.number_of_nodes(),
        "selected_source_side_road_2edge_edges": selected_road_graph.number_of_edges(),
        **terminal_minor_stats,
        "k_nearest_terminal_edges": K_NEAREST_TERMINAL_EDGES,
        **pruning_status,
        "selected_terminal_2edge_nodes": selected_terminal_graph.number_of_nodes(),
        "selected_terminal_2edge_edges": selected_terminal_graph.number_of_edges(),
        "transformers_total": int(len(transformers)),
        "transformers_in_selected_2edge": int(
            sum(selected_terminal_graph.nodes[n].get("transformer_count", 0) for n in selected_terminal_graph.nodes)
        ),
        "sources_total": int(len(substations)),
        "sources_in_selected_2edge": int(
            sum(selected_terminal_graph.nodes[n].get("source_count", 0) for n in selected_terminal_graph.nodes)
        ),
        "terminal_road_nodes": int(len(terminal_by_road_node)),
        "terminal_transformer_road_nodes": int(sum(1 for d in terminal_data.values() if d["kind"] in transformer_like)),
        "terminal_source_road_nodes": int(sum(1 for d in terminal_data.values() if d["kind"] in source_like)),
        "terminal_combined_source_transformer_road_nodes": int(
            sum(1 for d in terminal_data.values() if d["kind"] == "source_transformer")
        ),
        "raw_terminal_edges": int(len(raw_terminal_edges)),
        "corridor_transformer_transformer_edges": int(kind_counts["transformer-transformer"]),
        "corridor_edges_touching_source_terminal": source_touching_edges,
        "corridor_edges_between_transformer_only_terminals": transformer_only_edges,
        "corridor_total_length_m": float(raw_terminal_edges["length_m"].sum()) if len(raw_terminal_edges) else 0.0,
        "corridor_mean_length_m": float(raw_terminal_edges["length_m"].mean()) if len(raw_terminal_edges) else 0.0,
        "corridor_max_length_m": float(raw_terminal_edges["length_m"].max()) if len(raw_terminal_edges) else 0.0,
        "ilp_nodes": int(len(ilp_nodes)),
        "ilp_transformer_nodes": int(len(ilp_transformer_nodes)),
        "ilp_source_nodes": int(len(ilp_source_nodes)),
        "ilp_edges": int(len(ilp_edges)),
        "ilp_total_length_m": float(ilp_edges["length_m"].sum())
        if len(ilp_edges)
        else 0.0,
        "ilp_mean_length_m": float(ilp_edges["length_m"].mean())
        if len(ilp_edges)
        else 0.0,
        "ilp_max_length_m": float(ilp_edges["length_m"].max())
        if len(ilp_edges)
        else 0.0,
        "ilp_internal_transformers_on_edges": int(ilp_edges["edge_transformer_count"].sum())
        if len(ilp_edges)
        else 0,
        "ilp_internal_capacity_kva_on_edges": float(ilp_edges["edge_size_kva"].sum())
        if len(ilp_edges)
        else 0.0,
        "transformer_capacity_kva_in_selected_2edge": float(
            sum(selected_terminal_graph.nodes[n].get("size_kva", 0.0) for n in selected_terminal_graph.nodes)
        ),
    }
    SUMMARY_JSON_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    _write_summary(summary)
    return summary


def _write_summary(summary: dict[str, int | float]) -> None:
    lines = [
        "# P2U Terminal Corridor Preparation",
        "",
        "This is a temporary inspection artifact for the street-constrained optimization pipeline.",
        "Street junctions are used only to route corridor geometry. The corridor graph nodes are transformer/source terminals.",
        "",
        "## Outputs",
        "",
        f"- QGIS project: `{QGS_PATH.name}`",
        f"- GeoPackage: `{GPKG_PATH.name}`",
        "",
        "## Counts",
        "",
        f"- Full street graph nodes: `{summary['road_graph_nodes']}`",
        f"- Full street graph edges: `{summary['road_graph_edges']}`",
        f"- Source-side road 2-edge nodes kept before terminal graph: `{summary['selected_source_side_road_2edge_nodes']}`",
        f"- Source-side road 2-edge edges kept before terminal graph: `{summary['selected_source_side_road_2edge_edges']}`",
        f"- Terminal nodes before source-side 2-edge selection: `{summary['terminal_nodes_before_selection']}`",
        f"- Terminal edges before source-side 2-edge selection: `{summary['terminal_edges_before_selection']}`",
        f"- Direct road edges between terminals: `{summary['direct_terminal_road_edges']}`",
        f"- Non-terminal road components with at least two terminal boundaries: `{summary['nonterminal_components_with_2plus_boundary_terminals']}`",
        f"- Max boundary terminals on one non-terminal road component: `{summary['max_boundary_terminals_on_nonterminal_component']}`",
        f"- k-nearest pruning parameter: `{summary['k_nearest_terminal_edges']}`",
        f"- Pruned terminal nodes: `{summary['pruned_terminal_nodes']}`",
        f"- Pruned terminal edges: `{summary['pruned_terminal_edges']}`",
        f"- Entire pruned graph is source-side 2-edge-connected: `{summary['is_entire_pruned_graph_source_side_2edge']}`",
        f"- Source-side 2-edge nodes after pruning: `{summary['source_side_2edge_nodes_after_pruning']}`",
        f"- Source-side 2-edge edges after pruning: `{summary['source_side_2edge_edges_after_pruning']}`",
        f"- Selected source-side terminal 2-edge nodes: `{summary['selected_terminal_2edge_nodes']}`",
        f"- Selected source-side terminal 2-edge edges: `{summary['selected_terminal_2edge_edges']}`",
        f"- Transformers in selected component: `{summary['transformers_in_selected_2edge']}` / `{summary['transformers_total']}`",
        f"- Sources in selected component: `{summary['sources_in_selected_2edge']}` / `{summary['sources_total']}`",
        f"- Terminal road nodes: `{summary['terminal_road_nodes']}`",
        f"- Transformer terminal road nodes: `{summary['terminal_transformer_road_nodes']}`",
        f"- Source terminal road nodes: `{summary['terminal_source_road_nodes']}`",
        f"- Combined source-transformer terminal road nodes: `{summary['terminal_combined_source_transformer_road_nodes']}`",
        f"- Raw terminal-minor edges in selected component: `{summary['raw_terminal_edges']}`",
        f"- Transformer-transformer corridor edges: `{summary['corridor_transformer_transformer_edges']}`",
        f"- Corridor edges touching a source terminal: `{summary['corridor_edges_touching_source_terminal']}`",
        f"- Corridor edges between transformer-only terminals: `{summary['corridor_edges_between_transformer_only_terminals']}`",
        f"- Mean raw terminal-minor edge length: `{summary['corridor_mean_length_m']:.2f}` m",
        f"- Max raw terminal-minor edge length: `{summary['corridor_max_length_m']:.2f}` m",
        f"- ILP nodes after degree-2 terminal-chain contraction: `{summary['ilp_nodes']}`",
        f"- ILP transformer nodes: `{summary['ilp_transformer_nodes']}`",
        f"- ILP source nodes: `{summary['ilp_source_nodes']}`",
        f"- ILP edges: `{summary['ilp_edges']}`",
        f"- Mean ILP edge length: `{summary['ilp_mean_length_m']:.2f}` m",
        f"- Max ILP edge length: `{summary['ilp_max_length_m']:.2f}` m",
        f"- Internal transformers moved onto ILP edges: `{summary['ilp_internal_transformers_on_edges']}`",
        f"- Internal transformer capacity moved onto ILP edges: `{summary['ilp_internal_capacity_kva_on_edges']:.1f}` kVA",
        f"- Transformer capacity in selected component: `{summary['transformer_capacity_kva_in_selected_2edge']:.1f}` kVA",
        "",
        "## Method",
        "",
        "1. Attach transformers and HVMV sources to nearest street nodes.",
        "2. Remove terminal nodes from the road graph.",
        "3. For every connected component of non-terminal road nodes, connect all boundary terminals by shortest road-path distance through that component.",
        "4. Add direct terminal-terminal road edges where they exist.",
        f"5. Keep the symmetric `{K_NEAREST_TERMINAL_EDGES}` shortest incident terminal edges per terminal.",
        "6. Compute the source-side 2-edge-connected component on this pruned terminal graph after contracting source terminals.",
        "7. Contract degree-2 transformer/source-terminal chains into `ilp_edges`, preserving internal transformer count/capacity on the edge.",
        "",
        "The resulting graph is intended for visual validation before the ILP stage.",
    ]
    SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    summary = build_corridor_outputs()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
