from __future__ import annotations

import html
import json
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from shapely.geometry import LineString


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "power" / "better_grids" / "SFO" / "P2U"
OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo"
BACKBONE_SUMMARY = OUTPUT_DIR / "p2u_ilp_2edge_solution_summary.json"
BACKBONE_QGS = OUTPUT_DIR / "p2u_ilp_2edge_solution.qgs"
OUTPUT_GPKG = OUTPUT_DIR / "p2u_ilp_with_dijkstra_trees_3857.gpkg"
OUTPUT_QGS = OUTPUT_DIR / "p2u_ilp_with_dijkstra_trees.qgs"
SUMMARY_JSON = OUTPUT_DIR / "p2u_ilp_with_dijkstra_trees_summary.json"
SUMMARY_MD = OUTPUT_DIR / "p2u_ilp_with_dijkstra_trees_summary.md"

SOURCE_CRS = "EPSG:32610"
TARGET_CRS = "EPSG:3857"


def _clean(value) -> str:
    return str(value).strip()


def _read_layer(name: str) -> gpd.GeoDataFrame:
    gdf = gpd.read_file(DATA_DIR / f"{name}.shp")
    if gdf.crs is None:
        gdf = gdf.set_crs(SOURCE_CRS)
    return gdf


def _nearest_nodes(points: gpd.GeoDataFrame, street_nodes: gpd.GeoDataFrame) -> tuple[np.ndarray, np.ndarray]:
    point_geom = points.to_crs(street_nodes.crs).geometry.centroid
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
    for _, row in branches.iterrows():
        u = _clean(row.Node_A)
        v = _clean(row.Node_B)
        if u == v:
            continue
        length_m = float(row.geometry.length)
        if graph.has_edge(u, v):
            graph.edges[u, v]["length_m"] = min(graph.edges[u, v]["length_m"], length_m)
        else:
            graph.add_edge(u, v, length_m=length_m)
    return graph, node_pos


def _line_from_nodes(nodes: list[str], node_pos: dict[str, tuple[float, float]]) -> LineString:
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


def _terminal_road_nodes(backbone_gpkg: Path) -> set[str]:
    transformers = gpd.read_file(backbone_gpkg, layer="solution_transformer_nodes")
    sources = gpd.read_file(backbone_gpkg, layer="solution_source_nodes")
    road_nodes = set(transformers["road_node"].astype(str).map(_clean))
    source_road_col = "street_node" if "street_node" in sources.columns else "road_node"
    road_nodes.update(sources[source_road_col].astype(str).map(_clean))
    road_nodes.discard("__SOURCE__")
    return road_nodes


def build_tree_outputs() -> dict[str, int | float | str]:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if OUTPUT_GPKG.exists():
        try:
            OUTPUT_GPKG.unlink()
        except PermissionError:
            pass

    backbone_summary = json.loads(BACKBONE_SUMMARY.read_text(encoding="utf-8"))
    backbone_gpkg = Path(backbone_summary["output_gpkg"])
    backbone_edges = gpd.read_file(backbone_gpkg, layer="solution_edges")
    backbone_transformers = gpd.read_file(backbone_gpkg, layer="solution_transformer_nodes")
    backbone_sources = gpd.read_file(backbone_gpkg, layer="solution_source_nodes")

    branches = _read_layer("StreetMap_branches")
    street_nodes = _read_layer("StreetMap_nodes").copy()
    transformers = _read_layer("DistribTransf_N")
    street_nodes["Node"] = street_nodes["Node"].map(_clean)
    road_graph, node_pos = _build_road_graph(branches, street_nodes)

    backbone_road_nodes = _terminal_road_nodes(backbone_gpkg)
    missing_backbone_nodes = sorted(node for node in backbone_road_nodes if node not in road_graph)
    if missing_backbone_nodes:
        raise ValueError(f"{len(missing_backbone_nodes)} backbone road nodes are missing from the street graph")

    transformer_street_nodes, attachment_distances = _nearest_nodes(transformers, street_nodes)
    transformers = transformers.copy()
    transformers["transformer"] = transformers["Node"].map(_clean)
    transformers["street_node"] = list(map(_clean, transformer_street_nodes))
    transformers["nearest_street_distance_m"] = attachment_distances.astype(float)
    transformers["size_kva"] = transformers["Size_kVA"].astype(float)
    transformers["in_backbone"] = transformers["street_node"].isin(backbone_road_nodes)

    dist, paths = nx.multi_source_dijkstra(
        road_graph,
        sources=list(backbone_road_nodes),
        weight="length_m",
    )

    tree_edge_keys: set[tuple[str, str]] = set()
    attached_rows = []
    unattached_count = 0
    for _, row in transformers[~transformers["in_backbone"]].iterrows():
        street_node = row["street_node"]
        if street_node not in paths:
            unattached_count += 1
            continue
        path = paths[street_node]
        for u, v in zip(path[:-1], path[1:]):
            tree_edge_keys.add(tuple(sorted((u, v))))
        attached_rows.append(
            {
                "transformer": row["transformer"],
                "street_node": street_node,
                "nearest_backbone_distance_m": float(dist[street_node]),
                "size_kva": float(row["size_kva"]),
                "geometry": row.geometry,
            }
        )

    tree_edge_rows = []
    for u, v in sorted(tree_edge_keys):
        tree_edge_rows.append(
            {
                "a_node": u,
                "b_node": v,
                "length_m": float(road_graph.edges[u, v]["length_m"]),
                "geometry": _line_from_nodes([u, v], node_pos),
            }
        )
    tree_edges = gpd.GeoDataFrame(tree_edge_rows, geometry="geometry", crs=SOURCE_CRS)
    attached_transformers = gpd.GeoDataFrame(attached_rows, geometry="geometry", crs=transformers.crs)
    unattached_transformers = transformers[
        (~transformers["in_backbone"]) & (~transformers["street_node"].isin(paths.keys()))
    ].copy()

    backbone_edges.to_file(OUTPUT_GPKG, layer="backbone_edges", driver="GPKG")
    backbone_transformers.to_file(OUTPUT_GPKG, layer="backbone_transformer_nodes", driver="GPKG")
    backbone_sources.to_file(OUTPUT_GPKG, layer="backbone_source_nodes", driver="GPKG")
    if not tree_edges.empty:
        tree_edges.to_crs(TARGET_CRS).to_file(OUTPUT_GPKG, layer="dijkstra_tree_edges", driver="GPKG")
    if not attached_transformers.empty:
        attached_transformers.to_crs(TARGET_CRS).to_file(OUTPUT_GPKG, layer="tree_transformers", driver="GPKG")
    if not unattached_transformers.empty:
        unattached_transformers.to_crs(TARGET_CRS).to_file(OUTPUT_GPKG, layer="unattached_transformers", driver="GPKG")

    summary = {
        "backbone_solution_gpkg": str(backbone_gpkg),
        "backbone_edges": int(len(backbone_edges)),
        "backbone_transformer_nodes": int(len(backbone_transformers)),
        "source_nodes": int(len(backbone_sources)),
        "backbone_road_nodes": int(len(backbone_road_nodes)),
        "total_transformers": int(len(transformers)),
        "transformers_in_backbone_by_road_node": int(transformers["in_backbone"].sum()),
        "tree_transformers_attached": int(len(attached_transformers)),
        "tree_transformers_unattached": int(unattached_count),
        "dijkstra_tree_edges": int(len(tree_edges)),
        "dijkstra_tree_length_m": float(tree_edges["length_m"].sum()) if not tree_edges.empty else 0.0,
        "tree_transformer_capacity_kva": float(attached_transformers["size_kva"].sum()) if not attached_transformers.empty else 0.0,
        "mean_distance_to_backbone_m": float(attached_transformers["nearest_backbone_distance_m"].mean())
        if not attached_transformers.empty
        else 0.0,
        "max_distance_to_backbone_m": float(attached_transformers["nearest_backbone_distance_m"].max())
        if not attached_transformers.empty
        else 0.0,
        "output_gpkg": str(OUTPUT_GPKG),
    }
    SUMMARY_JSON.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    SUMMARY_MD.write_text(_summary_md(summary), encoding="utf-8")
    write_qgs()
    return summary


def _summary_md(summary: dict) -> str:
    return "\n".join(
        [
            "# P2U ILP Backbone With Dijkstra Trees",
            "",
            f"- Backbone edges: `{summary['backbone_edges']}`",
            f"- Backbone transformer nodes: `{summary['backbone_transformer_nodes']}`",
            f"- Source nodes: `{summary['source_nodes']}`",
            f"- Total transformers: `{summary['total_transformers']}`",
            f"- Transformers in backbone by road node: `{summary['transformers_in_backbone_by_road_node']}`",
            f"- Tree transformers attached: `{summary['tree_transformers_attached']}`",
            f"- Tree transformers unattached: `{summary['tree_transformers_unattached']}`",
            f"- Dijkstra tree edges: `{summary['dijkstra_tree_edges']}`",
            f"- Dijkstra tree length: `{summary['dijkstra_tree_length_m']:.3f}` m",
            f"- Tree transformer capacity: `{summary['tree_transformer_capacity_kva']:.1f}` kVA",
            f"- Mean/max distance to backbone: `{summary['mean_distance_to_backbone_m']:.3f}` m / `{summary['max_distance_to_backbone_m']:.3f}` m",
            "",
            "The tree layer is the union of shortest road paths from each non-backbone transformer to the optimized 2-edge backbone.",
        ]
    ) + "\n"


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


def _map_layer(layer: str, name: str, layer_id: str, geometry: str, renderer: str, checked: str = "1") -> tuple[str, str]:
    tree = f'      <layer-tree-layer checked="{checked}" id="{layer_id}" name="{html.escape(name)}" source="./{OUTPUT_GPKG.name}|layername={layer}"/>'
    project = f"""  <maplayer type="vector" geometry="{geometry}">
    <id>{layer_id}</id>
    <datasource>./{OUTPUT_GPKG.name}|layername={layer}</datasource>
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
    return tree, project


def write_qgs() -> None:
    layer_specs = [
        ("backbone_edges", "ILP 2-edge backbone", "backbone_edges_0de91cb1", "Line", _line_symbol("174,56,132,255", "0.65"), "1"),
        ("dijkstra_tree_edges", "Dijkstra 1-comp trees", "dijkstra_tree_edges_78fd5ba3", "Line", _line_symbol("64,118,184,210", "0.35"), "1"),
        ("backbone_transformer_nodes", "Backbone transformers", "backbone_transformer_nodes_b6b97241", "Point", _point_symbol("255,255,255,0", "174,56,132,255", "1.2", False), "1"),
        ("tree_transformers", "Tree transformers", "tree_transformers_feb29bc4", "Point", _point_symbol("255,255,255,0", "64,118,184,255", "0.9", False), "1"),
        ("backbone_source_nodes", "Sources", "backbone_source_nodes_72d34642", "Point", _point_symbol("18,18,18,255", "255,255,255,255", "2.6", True), "1"),
    ]
    trees = []
    projects = []
    for spec in layer_specs:
        try:
            gpd.read_file(OUTPUT_GPKG, layer=spec[0])
        except Exception:
            continue
        tree, project = _map_layer(*spec)
        trees.append(tree)
        projects.append(project)

    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="P2U ILP backbone with Dijkstra trees">
  <homePath path="."/>
  <title>P2U ILP backbone with Dijkstra trees</title>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="P2U optimized network">
{chr(10).join(trees)}
    </layer-tree-group>
  </layer-tree-group>
  <projectlayers>
{chr(10).join(projects)}
  </projectlayers>
</qgis>
"""
    OUTPUT_QGS.write_text(qgs, encoding="utf-8")


def main() -> None:
    summary = build_tree_outputs()
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
