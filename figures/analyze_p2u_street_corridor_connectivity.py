from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "power" / "better_grids" / "SFO" / "P2U"
OUTPUT_DIR = ROOT / "outputs" / "p2u_city_corridor_feasibility"
SOURCE_CRS = "EPSG:32610"


def _read_layer(name: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(DATA_DIR / f"{name}.shp")
    if gdf.crs is None:
        gdf = gdf.set_crs(SOURCE_CRS)
    return gdf


def _clean(value) -> str:
    return str(value).strip()


def _nearest_nodes(points: gpd.GeoDataFrame, street_nodes: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    points = points.to_crs(street_nodes.crs)
    point_geom = points.geometry.centroid
    node_geom = street_nodes.geometry
    tree = cKDTree(np.column_stack([node_geom.x.to_numpy(), node_geom.y.to_numpy()]))
    distances, indices = tree.query(np.column_stack([point_geom.x.to_numpy(), point_geom.y.to_numpy()]), k=1)
    return street_nodes.iloc[indices]["Node"].to_numpy(), distances


def _build_street_graph(branches: gpd.GeoDataFrame) -> nx.Graph:
    graph = nx.Graph()
    duplicate_edges = 0
    self_loops = 0
    for _, row in branches.iterrows():
        u = _clean(row.Node_A)
        v = _clean(row.Node_B)
        if u == v:
            self_loops += 1
            continue
        length_m = float(row.geometry.length)
        if graph.has_edge(u, v):
            duplicate_edges += 1
            graph.edges[u, v]["duplicate_count"] += 1
            graph.edges[u, v]["length_m"] = min(graph.edges[u, v]["length_m"], length_m)
        else:
            graph.add_edge(u, v, length_m=length_m, duplicate_count=1)
    graph.graph["duplicate_edges"] = duplicate_edges
    graph.graph["self_loops"] = self_loops
    return graph


def _source_side_2edge_nodes(graph: nx.Graph, source_nodes: set[str]) -> tuple[set[str], dict[str, int], int]:
    """Return nodes edge-2-connected to at least one source-side street node."""
    source_components = []
    comp_lookup = {}
    for comp_id, comp_nodes in enumerate(nx.connected_components(graph)):
        comp_nodes = set(comp_nodes)
        if comp_nodes & source_nodes:
            source_components.append(comp_nodes)
            for node in comp_nodes:
                comp_lookup[node] = comp_id

    two_edge_nodes = set()
    two_edge_component_sizes = {}
    source_connected_component_count = len(source_components)
    for comp_nodes in source_components:
        sub = graph.subgraph(comp_nodes).copy()
        for block_id, block_nodes in enumerate(nx.k_edge_components(sub, k=2)):
            block_nodes = set(block_nodes)
            if block_nodes & source_nodes:
                two_edge_nodes.update(block_nodes)
                two_edge_component_sizes[f"source_block_{len(two_edge_component_sizes)}"] = len(block_nodes)
    return two_edge_nodes, comp_lookup, source_connected_component_count


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    branches = _read_layer("StreetMap_branches")
    street_nodes = _read_layer("StreetMap_nodes")
    transformers = _read_layer("DistribTransf_N")
    substations = _read_layer("HVMVSubstation_N")

    street_nodes = street_nodes.copy()
    street_nodes["Node"] = street_nodes["Node"].map(_clean)

    graph = _build_street_graph(branches)

    transformer_street_nodes, transformer_distances = _nearest_nodes(transformers, street_nodes)
    source_street_nodes, source_distances = _nearest_nodes(substations, street_nodes)
    source_street_node_set = set(map(_clean, source_street_nodes))

    two_edge_nodes, comp_lookup, source_connected_component_count = _source_side_2edge_nodes(
        graph, source_street_node_set
    )

    transformer_records = []
    counts = Counter()
    capacity_by_status = Counter()
    source_component_transformer_counts = Counter()
    component_transformer_counts = Counter()
    distance_by_status = defaultdict(list)

    for i, (_, row) in enumerate(transformers.iterrows()):
        street_node = _clean(transformer_street_nodes[i])
        connected_to_source = street_node in comp_lookup
        two_edge_to_source = street_node in two_edge_nodes
        if not connected_to_source:
            status = "not_source_connected"
        elif two_edge_to_source:
            status = "street_2edge_to_source"
        else:
            status = "street_1edge_to_source"

        counts[status] += 1
        size_kva = float(row.get("Size_kVA", 0.0) or 0.0)
        capacity_by_status[status] += size_kva
        if connected_to_source:
            source_component_transformer_counts[comp_lookup[street_node]] += 1
        component_transformer_counts[street_node] += 1
        distance_by_status[status].append(float(transformer_distances[i]))

        transformer_records.append(
            {
                "transformer": _clean(row.get("Code", row.get("Node", i))),
                "original_node": _clean(row.get("Node", "")),
                "size_kva": size_kva,
                "street_node": street_node,
                "nearest_street_distance_m": float(transformer_distances[i]),
                "status": status,
                "connected_to_source_street_component": connected_to_source,
                "edge_2_connected_to_source_on_street_graph": two_edge_to_source,
            }
        )

    connected_components = list(nx.connected_components(graph))
    bridge_count = len(list(nx.bridges(graph)))
    k2_components = list(nx.k_edge_components(graph, k=2))
    source_k2_components = [
        len(set(block) & source_street_node_set)
        for block in k2_components
        if set(block) & source_street_node_set
    ]

    summary = {
        "definition": (
            "A transformer is street_2edge_to_source if its nearest StreetMap node is in the "
            "same 2-edge-connected block as at least one nearest HVMVSubstation street node. "
            "A street_1edge_to_source transformer is source-connected through the city corridor graph, "
            "but every possible street-corridor route to a source crosses at least one bridge."
        ),
        "street_graph": {
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "raw_branch_features": len(branches),
            "duplicate_edge_features_collapsed": graph.graph["duplicate_edges"],
            "self_loop_features_skipped": graph.graph["self_loops"],
            "connected_components": len(connected_components),
            "largest_component_nodes": max(len(c) for c in connected_components),
            "bridges": bridge_count,
            "k2_components": len(k2_components),
            "source_connected_components": source_connected_component_count,
            "source_street_nodes": len(source_street_node_set),
            "source_k2_blocks": len(source_k2_components),
        },
        "transformers": {
            "total": len(transformers),
            "street_2edge_to_source": counts["street_2edge_to_source"],
            "street_1edge_to_source": counts["street_1edge_to_source"],
            "not_source_connected": counts["not_source_connected"],
            "cannot_be_made_2edge_using_only_street_corridors": (
                counts["street_1edge_to_source"] + counts["not_source_connected"]
            ),
        },
        "transformer_capacity_kva": {
            "total": float(sum(capacity_by_status.values())),
            "street_2edge_to_source": float(capacity_by_status["street_2edge_to_source"]),
            "street_1edge_to_source": float(capacity_by_status["street_1edge_to_source"]),
            "not_source_connected": float(capacity_by_status["not_source_connected"]),
            "cannot_be_made_2edge_using_only_street_corridors": float(
                capacity_by_status["street_1edge_to_source"]
                + capacity_by_status["not_source_connected"]
            ),
        },
        "nearest_attachment_distance_m": {
            status: {
                "count": len(values),
                "mean": float(np.mean(values)) if values else 0.0,
                "median": float(np.median(values)) if values else 0.0,
                "max": float(np.max(values)) if values else 0.0,
            }
            for status, values in distance_by_status.items()
        },
        "source_attachment_distance_m": {
            "count": len(source_distances),
            "mean": float(np.mean(source_distances)) if len(source_distances) else 0.0,
            "median": float(np.median(source_distances)) if len(source_distances) else 0.0,
            "max": float(np.max(source_distances)) if len(source_distances) else 0.0,
        },
    }

    (OUTPUT_DIR / "p2u_street_corridor_connectivity_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    gpd.GeoDataFrame(transformer_records, geometry=transformers.geometry, crs=transformers.crs).to_file(
        OUTPUT_DIR / "p2u_transformer_street_corridor_status.gpkg",
        layer="transformer_street_corridor_status",
        driver="GPKG",
    )

    md_lines = [
        "# P2U Street-Corridor 2-Edge Feasibility",
        "",
        "## Definition",
        "",
        summary["definition"],
        "",
        "## Street Graph",
        "",
        f"- Nodes: `{summary['street_graph']['nodes']}`",
        f"- Edges after duplicate collapse: `{summary['street_graph']['edges']}`",
        f"- Raw branch features: `{summary['street_graph']['raw_branch_features']}`",
        f"- Duplicate edge features collapsed: `{summary['street_graph']['duplicate_edge_features_collapsed']}`",
        f"- Connected components: `{summary['street_graph']['connected_components']}`",
        f"- Bridges: `{summary['street_graph']['bridges']}`",
        f"- 2-edge-connected blocks: `{summary['street_graph']['k2_components']}`",
        f"- Source-side connected components: `{summary['street_graph']['source_connected_components']}`",
        "",
        "## Transformers",
        "",
        f"- Total transformers: `{summary['transformers']['total']}`",
        f"- 2-edge-connected to a source through street corridors: `{summary['transformers']['street_2edge_to_source']}`",
        f"- Source-connected but only 1-edge-connected through street corridors: `{summary['transformers']['street_1edge_to_source']}`",
        f"- Not connected to any source-side street component: `{summary['transformers']['not_source_connected']}`",
        f"- Cannot be made 2-edge-connected using only existing street corridors: `{summary['transformers']['cannot_be_made_2edge_using_only_street_corridors']}`",
        f"- Total transformer capacity: `{summary['transformer_capacity_kva']['total']:.1f}` kVA",
        f"- Capacity 2-edge-connected through street corridors: `{summary['transformer_capacity_kva']['street_2edge_to_source']:.1f}` kVA",
        f"- Capacity source-connected but only 1-edge-connected: `{summary['transformer_capacity_kva']['street_1edge_to_source']:.1f}` kVA",
        f"- Capacity impossible to make 2-edge-connected using only existing street corridors: `{summary['transformer_capacity_kva']['cannot_be_made_2edge_using_only_street_corridors']:.1f}` kVA",
        "",
        "## Attachment Distances",
        "",
    ]
    for status, stats in sorted(summary["nearest_attachment_distance_m"].items()):
        md_lines.append(
            f"- {status}: count `{stats['count']}`, mean `{stats['mean']:.2f}` m, "
            f"median `{stats['median']:.2f}` m, max `{stats['max']:.2f}` m"
        )
    md_lines.extend(
        [
            "",
            "## Notes",
            "",
            "- This is a geometric/corridor feasibility analysis, not an electrical power-flow feasibility test.",
            "- Transformers and substations are attached to their nearest `StreetMap_nodes` point.",
            "- The street graph is treated as a simple graph; duplicate street branch features are collapsed.",
        ]
    )
    (OUTPUT_DIR / "p2u_street_corridor_connectivity_summary.md").write_text(
        "\n".join(md_lines) + "\n", encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
