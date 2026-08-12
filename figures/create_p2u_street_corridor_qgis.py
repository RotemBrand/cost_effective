from __future__ import annotations

import html
import json
from collections import Counter
from pathlib import Path

import geopandas as gpd
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "power" / "better_grids" / "SFO" / "P2U"
QGIS_DIR = ROOT / "outputs" / "qgis" / "P2U"
GPKG_PATH = QGIS_DIR / "P2U_street_corridor_connectivity_3857.gpkg"
QGS_PATH = QGIS_DIR / "P2U_street_corridor_connectivity.qgs"
README_PATH = QGIS_DIR / "README_STREET_CORRIDOR_CONNECTIVITY.md"
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
    points = points.to_crs(street_nodes.crs)
    point_geom = points.geometry.centroid
    node_geom = street_nodes.geometry
    tree = cKDTree(np.column_stack([node_geom.x.to_numpy(), node_geom.y.to_numpy()]))
    distances, indices = tree.query(np.column_stack([point_geom.x.to_numpy(), point_geom.y.to_numpy()]), k=1)
    return street_nodes.iloc[indices]["Node"].to_numpy(), distances


def _build_street_graph(branches: gpd.GeoDataFrame) -> nx.Graph:
    graph = nx.Graph()
    for _, row in branches.iterrows():
        u = _clean(row.Node_A)
        v = _clean(row.Node_B)
        if u == v:
            continue
        graph.add_edge(u, v, length_m=float(row.geometry.length))
    return graph


def _source_connectivity_sets(
    graph: nx.Graph, source_street_nodes: set[str]
) -> tuple[set[str], set[str], dict[str, int]]:
    source_connected_nodes = set()
    source_component_id_by_node: dict[str, int] = {}
    for component_id, component in enumerate(nx.connected_components(graph)):
        component = set(component)
        if component & source_street_nodes:
            source_connected_nodes.update(component)
            for node in component:
                source_component_id_by_node[node] = component_id

    source_side_2edge_nodes = set()
    for component in nx.connected_components(graph):
        component = set(component)
        if not component & source_street_nodes:
            continue
        subgraph = graph.subgraph(component).copy()
        for block in nx.k_edge_components(subgraph, k=2):
            block = set(block)
            if block & source_street_nodes:
                source_side_2edge_nodes.update(block)

    return source_connected_nodes, source_side_2edge_nodes, source_component_id_by_node


def _write_gpkg_layer(gdf: gpd.GeoDataFrame, layer: str) -> None:
    if gdf.empty:
        return
    gdf.to_crs(TARGET_CRS).to_file(GPKG_PATH, layer=layer, driver="GPKG")


def build_layers() -> dict[str, int | float]:
    QGIS_DIR.mkdir(parents=True, exist_ok=True)
    if GPKG_PATH.exists():
        GPKG_PATH.unlink()

    branches = _read_layer("StreetMap_branches")
    street_nodes = _read_layer("StreetMap_nodes").copy()
    transformers = _read_layer("DistribTransf_N")
    substations = _read_layer("HVMVSubstation_N")

    street_nodes["Node"] = street_nodes["Node"].map(_clean)
    graph = _build_street_graph(branches)
    transformer_street_nodes, transformer_distances = _nearest_nodes(transformers, street_nodes)
    source_street_nodes, source_distances = _nearest_nodes(substations, street_nodes)
    source_street_node_set = set(map(_clean, source_street_nodes))

    source_connected_nodes, source_side_2edge_nodes, source_component_id_by_node = _source_connectivity_sets(
        graph, source_street_node_set
    )

    street_rows = branches.copy()
    classes = []
    for _, row in street_rows.iterrows():
        u = _clean(row.Node_A)
        v = _clean(row.Node_B)
        if u in source_side_2edge_nodes and v in source_side_2edge_nodes:
            classes.append("source_side_2edge_corridor")
        elif u in source_connected_nodes and v in source_connected_nodes:
            classes.append("source_connected_1edge_corridor")
        else:
            classes.append("not_source_connected_corridor")
    street_rows["corridor_status"] = classes
    street_rows["a_node"] = street_rows["Node_A"].map(_clean)
    street_rows["b_node"] = street_rows["Node_B"].map(_clean)
    street_rows["length_m"] = street_rows.geometry.length.astype(float)
    street_rows = street_rows[["corridor_status", "a_node", "b_node", "length_m", "geometry"]]

    transformer_rows = transformers.copy()
    transformer_status = []
    transformer_street_node_clean = []
    transformer_component_ids = []
    for street_node in map(_clean, transformer_street_nodes):
        transformer_street_node_clean.append(street_node)
        transformer_component_ids.append(source_component_id_by_node.get(street_node, -1))
        if street_node in source_side_2edge_nodes:
            transformer_status.append("street_2edge_to_source")
        elif street_node in source_connected_nodes:
            transformer_status.append("street_1edge_to_source")
        else:
            transformer_status.append("not_source_connected")
    transformer_rows["transformer"] = transformer_rows["Node"].map(_clean)
    transformer_rows["street_node"] = transformer_street_node_clean
    transformer_rows["source_component_id"] = transformer_component_ids
    transformer_rows["nearest_street_distance_m"] = transformer_distances.astype(float)
    transformer_rows["status"] = transformer_status
    transformer_rows["size_kva"] = transformer_rows["Size_kVA"].astype(float)
    transformer_rows = transformer_rows[
        [
            "transformer",
            "street_node",
            "source_component_id",
            "nearest_street_distance_m",
            "status",
            "size_kva",
            "NomV_kV",
            "Phases",
            "geometry",
        ]
    ]

    source_rows = substations.copy()
    source_rows["source"] = source_rows["Node"].map(_clean)
    source_rows["street_node"] = list(map(_clean, source_street_nodes))
    source_rows["nearest_street_distance_m"] = source_distances.astype(float)
    source_rows = source_rows[["source", "street_node", "nearest_street_distance_m", "geometry"]]

    _write_gpkg_layer(
        street_rows[street_rows["corridor_status"] == "source_side_2edge_corridor"],
        "street_2edge_corridors",
    )
    _write_gpkg_layer(
        street_rows[street_rows["corridor_status"] == "source_connected_1edge_corridor"],
        "street_1edge_corridors",
    )
    _write_gpkg_layer(
        street_rows[street_rows["corridor_status"] == "not_source_connected_corridor"],
        "street_not_source_connected_corridors",
    )
    _write_gpkg_layer(
        transformer_rows[transformer_rows["status"] == "street_2edge_to_source"],
        "transformers_2edge",
    )
    _write_gpkg_layer(
        transformer_rows[transformer_rows["status"] == "street_1edge_to_source"],
        "transformers_1edge",
    )
    _write_gpkg_layer(
        transformer_rows[transformer_rows["status"] == "not_source_connected"],
        "transformers_not_source_connected",
    )
    _write_gpkg_layer(source_rows, "sources")

    counts = Counter(classes)
    transformer_counts = Counter(transformer_status)
    summary = {
        "street_nodes": graph.number_of_nodes(),
        "street_edges": graph.number_of_edges(),
        "street_2edge_corridors": counts["source_side_2edge_corridor"],
        "street_1edge_corridors": counts["source_connected_1edge_corridor"],
        "street_not_source_connected_corridors": counts["not_source_connected_corridor"],
        "transformers_total": len(transformer_rows),
        "transformers_2edge": transformer_counts["street_2edge_to_source"],
        "transformers_1edge": transformer_counts["street_1edge_to_source"],
        "transformers_not_source_connected": transformer_counts["not_source_connected"],
        "transformer_capacity_kva_total": float(transformer_rows["size_kva"].sum()),
        "transformer_capacity_kva_2edge": float(
            transformer_rows.loc[transformer_rows["status"] == "street_2edge_to_source", "size_kva"].sum()
        ),
        "transformer_capacity_kva_1edge": float(
            transformer_rows.loc[transformer_rows["status"] == "street_1edge_to_source", "size_kva"].sum()
        ),
    }
    return summary


def _layer_tree(ids: dict[str, str]) -> str:
    layers = [
        ("Street 2-edge corridors", ids["street_2edge_corridors"], "1"),
        ("Street 1-edge corridors", ids["street_1edge_corridors"], "1"),
        ("Street off-source corridors", ids["street_not_source_connected_corridors"], "0"),
        ("Transformers 2-edge", ids["transformers_2edge"], "1"),
        ("Transformers 1-edge", ids["transformers_1edge"], "1"),
        ("Sources", ids["sources"], "1"),
    ]
    entries = "\n".join(
        f'      <layer-tree-layer checked="{checked}" id="{layer_id}" name="{html.escape(name)}" source="./{GPKG_PATH.name}|layername={layer}"/>'
        for (name, layer_id, checked), layer in zip(layers, ids)
    )
    return f"""  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="Street corridor connectivity">
{entries}
    </layer-tree-group>
  </layer-tree-group>"""


def _line_symbol(color: str, width: str, dash: str = "solid") -> str:
    return f"""      <renderer-v2 type="singleSymbol" enableorderby="0" forceraster="0" referencescale="-1" symbollevels="0">
        <symbols>
          <symbol alpha="1" type="line" name="0" force_rhr="0" clip_to_extent="1" is_animated="0" frame_rate="10">
            <layer enabled="1" id="{{2f1cb3e4-2c04-4c6f-9d5d-515c6e8fb273}}" locked="0" pass="0" class="SimpleLine">
              <Option type="Map">
                <Option name="line_color" type="QString" value="{color}"/>
                <Option name="line_style" type="QString" value="{dash}"/>
                <Option name="line_width" type="QString" value="{width}"/>
                <Option name="line_width_unit" type="QString" value="MM"/>
              </Option>
            </layer>
          </symbol>
        </symbols>
      </renderer-v2>"""


def _point_symbol(color: str, outline: str, size: str, fill: bool = True) -> str:
    fill_color = color if fill else "255,255,255,0"
    return f"""      <renderer-v2 type="singleSymbol" enableorderby="0" forceraster="0" referencescale="-1" symbollevels="0">
        <symbols>
          <symbol alpha="1" type="marker" name="0" force_rhr="0" clip_to_extent="1" is_animated="0" frame_rate="10">
            <layer enabled="1" id="{{0f6fdc49-d520-4669-bf7e-48e9e35899d8}}" locked="0" pass="0" class="SimpleMarker">
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
    return f"""  <maplayer type="vector" styleCategories="AllStyleCategories" refreshOnNotifyEnabled="0" autoRefreshEnabled="0" autoRefreshTime="0" minScale="100000000" maxScale="0" simplifyMaxScale="1" hasScaleBasedVisibilityFlag="0" readOnly="0" geometry="{geometry}">
    <id>{layer_id}</id>
    <datasource>./{GPKG_PATH.name}|layername={layer}</datasource>
    <layername>{html.escape(name)}</layername>
    <srs>
      <spatialrefsys nativeFormat="Wkt">
        <wkt>PROJCRS[&quot;WGS 84 / Pseudo-Mercator&quot;,BASEGEOGCRS[&quot;WGS 84&quot;,ENSEMBLE[&quot;World Geodetic System 1984 ensemble&quot;,MEMBER[&quot;World Geodetic System 1984 (Transit)&quot;],MEMBER[&quot;World Geodetic System 1984 (G730)&quot;],MEMBER[&quot;World Geodetic System 1984 (G873)&quot;],MEMBER[&quot;World Geodetic System 1984 (G1150)&quot;],MEMBER[&quot;World Geodetic System 1984 (G1674)&quot;],MEMBER[&quot;World Geodetic System 1984 (G1762)&quot;],MEMBER[&quot;World Geodetic System 1984 (G2139)&quot;],ELLIPSOID[&quot;WGS 84&quot;,6378137,298.257223563,LENGTHUNIT[&quot;metre&quot;,1]],ENSEMBLEACCURACY[2.0]],PRIMEM[&quot;Greenwich&quot;,0,ANGLEUNIT[&quot;degree&quot;,0.0174532925199433]],ID[&quot;EPSG&quot;,4326]],CONVERSION[&quot;Popular Visualisation Pseudo-Mercator&quot;,METHOD[&quot;Popular Visualisation Pseudo Mercator&quot;,ID[&quot;EPSG&quot;,1024]],PARAMETER[&quot;Latitude of natural origin&quot;,0,ANGLEUNIT[&quot;degree&quot;,0.0174532925199433],ID[&quot;EPSG&quot;,8801]],PARAMETER[&quot;Longitude of natural origin&quot;,0,ANGLEUNIT[&quot;degree&quot;,0.0174532925199433],ID[&quot;EPSG&quot;,8802]],PARAMETER[&quot;False easting&quot;,0,LENGTHUNIT[&quot;metre&quot;,1],ID[&quot;EPSG&quot;,8806]],PARAMETER[&quot;False northing&quot;,0,LENGTHUNIT[&quot;metre&quot;,1],ID[&quot;EPSG&quot;,8807]]],CS[Cartesian,2],AXIS[&quot;easting (X)&quot;,east,ORDER[1],LENGTHUNIT[&quot;metre&quot;,1]],AXIS[&quot;northing (Y)&quot;,north,ORDER[2],LENGTHUNIT[&quot;metre&quot;,1]],USAGE[SCOPE[&quot;Web mapping and visualisation.&quot;],AREA[&quot;World between 85.06°S and 85.06°N.&quot;],BBOX[-85.06,-180,85.06,180]],ID[&quot;EPSG&quot;,3857]]</wkt>
        <proj4>+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs</proj4>
        <srsid>3857</srsid>
        <srid>3857</srid>
        <authid>EPSG:3857</authid>
        <description>WGS 84 / Pseudo-Mercator</description>
        <projectionacronym>merc</projectionacronym>
        <ellipsoidacronym>EPSG:7030</ellipsoidacronym>
        <geographicflag>false</geographicflag>
      </spatialrefsys>
    </srs>
    <provider encoding="UTF-8">ogr</provider>
{renderer}
    <labeling type="simple"/>
  </maplayer>"""


def write_qgs() -> None:
    ids = {
        "street_2edge_corridors": "street_2edge_corridors_42f780a1",
        "street_1edge_corridors": "street_1edge_corridors_744bfbb2",
        "street_not_source_connected_corridors": "street_not_source_connected_corridors_72e9bde8",
        "transformers_2edge": "transformers_2edge_e0f8b635",
        "transformers_1edge": "transformers_1edge_c351159d",
        "sources": "sources_c6735728",
    }
    layers = [
        _map_layer(
            "street_2edge_corridors",
            "Street 2-edge corridors",
            ids["street_2edge_corridors"],
            "Line",
            _line_symbol("36,128,86,255", "0.28"),
        ),
        _map_layer(
            "street_1edge_corridors",
            "Street 1-edge corridors",
            ids["street_1edge_corridors"],
            "Line",
            _line_symbol("203,64,60,255", "0.42"),
        ),
        _map_layer(
            "street_not_source_connected_corridors",
            "Street off-source corridors",
            ids["street_not_source_connected_corridors"],
            "Line",
            _line_symbol("145,145,145,170", "0.20", "dash"),
        ),
        _map_layer(
            "transformers_2edge",
            "Transformers 2-edge",
            ids["transformers_2edge"],
            "Point",
            _point_symbol("255,255,255,0", "36,128,86,255", "1.1", fill=False),
        ),
        _map_layer(
            "transformers_1edge",
            "Transformers 1-edge",
            ids["transformers_1edge"],
            "Point",
            _point_symbol("255,255,255,0", "203,64,60,255", "1.35", fill=False),
        ),
        _map_layer(
            "sources",
            "Sources",
            ids["sources"],
            "Point",
            _point_symbol("18,18,18,255", "255,255,255,255", "2.4", fill=True),
        ),
    ]
    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="P2U street corridor connectivity">
  <homePath path="."/>
  <title>P2U street corridor connectivity</title>
{_layer_tree(ids)}
  <projectlayers>
{chr(10).join(layers)}
  </projectlayers>
</qgis>
"""
    QGS_PATH.write_text(qgs, encoding="utf-8")


def write_readme(summary: dict[str, int | float]) -> None:
    README_PATH.write_text(
        "\n".join(
            [
                "# P2U Street Corridor Connectivity QGIS",
                "",
                "Generated files:",
                "",
                f"- `{GPKG_PATH.name}`",
                f"- `{QGS_PATH.name}`",
                "",
                "Layers:",
                "",
                "- `street_2edge_corridors`: street corridors in source-side 2-edge-connected blocks.",
                "- `street_1edge_corridors`: street corridors source-connected only through at least one bridge.",
                "- `street_not_source_connected_corridors`: street corridors outside the source-connected street component.",
                "- `transformers_2edge`: transformers attached to source-side 2-edge street nodes.",
                "- `transformers_1edge`: transformers attached to source-connected but 1-edge street nodes.",
                "- `sources`: HVMV substations attached to street nodes.",
                "",
                "Summary:",
                "",
                f"- Street edges in source-side 2-edge corridors: `{summary['street_2edge_corridors']}`",
                f"- Street edges in source-connected 1-edge corridors: `{summary['street_1edge_corridors']}`",
                f"- Street edges not source-connected: `{summary['street_not_source_connected_corridors']}`",
                f"- Transformers in source-side 2-edge corridors: `{summary['transformers_2edge']}`",
                f"- Transformers in source-connected 1-edge corridors: `{summary['transformers_1edge']}`",
                f"- Transformer capacity in source-side 2-edge corridors: `{summary['transformer_capacity_kva_2edge']:.1f}` kVA",
                f"- Transformer capacity in source-connected 1-edge corridors: `{summary['transformer_capacity_kva_1edge']:.1f}` kVA",
                "",
                "This is a corridor-topology feasibility view. It does not validate electrical constraints.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    summary = build_layers()
    write_qgs()
    write_readme(summary)
    (QGIS_DIR / "P2U_street_corridor_connectivity_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    print(f"wrote {GPKG_PATH.name}")
    print(f"wrote {QGS_PATH.name}")


if __name__ == "__main__":
    main()
