from __future__ import annotations

import html
import re
import sys
from pathlib import Path

import geopandas as gpd
import networkx as nx
import pyogrio

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.p2u_full_mv_reliability import (
    SOURCE_CRS,
    _clean_node,
    _edge_key,
    build_full_mv_graph,
)


QGIS_DIR = ROOT / "outputs" / "qgis" / "P2U"
GPKG_PATH = QGIS_DIR / "P2U_MV_physical_network_3857.gpkg"
QGS_PATH = QGIS_DIR / "P2U_MV_physical_network_open_ties.qgs"
README_PATH = QGIS_DIR / "README_MV_PHYSICAL_OPEN_TIES.md"

SOURCE_NODE = "__SOURCE__"
BRIDGE_LAYER = "mv_bridge_lines"
COMPONENT_LAYER = "mv_2edge_component_lines"
BRIDGE_LAYER_ID = "MV_bridge_lines_8f3c6e21"
COMPONENT_LAYER_ID = "MV_2edge_component_lines_d4e98712"
TIE_LAYER_ID = "MV_normally_open_ties_ad74edf3"


def main() -> None:
    graph, mv_lines, _ = build_full_mv_graph()
    bridge_gdf, component_gdf = build_connectivity_layers(graph, mv_lines)
    write_layers(bridge_gdf, component_gdf)
    update_qgs_project(component_ids=sorted(component_gdf["component_id"].unique().tolist()))
    update_readme(len(bridge_gdf), len(component_gdf), component_gdf["component_id"].nunique())
    print(
        f"wrote {BRIDGE_LAYER} ({len(bridge_gdf)} features) and "
        f"{COMPONENT_LAYER} ({len(component_gdf)} features, "
        f"{component_gdf['component_id'].nunique()} components)"
    )


def build_connectivity_layers(
    graph: nx.Graph,
    mv_lines: gpd.GeoDataFrame,
) -> tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]:
    """Classify full MV physical lines after contracting all source nodes."""
    sources = set(graph.graph.get("sources", []))
    contracted = nx.relabel_nodes(graph, {source: SOURCE_NODE for source in sources}, copy=True)
    contracted.remove_edges_from(nx.selfloop_edges(contracted))

    bridges = {_edge_key(u, v) for u, v in nx.bridges(contracted)}
    without_bridges = contracted.copy()
    without_bridges.remove_edges_from(bridges)
    components = [set(comp) for comp in nx.connected_components(without_bridges)]
    node_to_component = {
        node: component_id
        for component_id, comp in enumerate(components)
        for node in comp
    }
    component_sizes = {component_id: len(comp) for component_id, comp in enumerate(components)}
    component_edge_counts = {component_id: 0 for component_id in range(len(components))}
    for u, v in without_bridges.edges():
        component_id = node_to_component[u]
        if component_id == node_to_component[v]:
            component_edge_counts[component_id] += 1

    rows = mv_lines.copy()
    if rows.crs is None:
        rows = rows.set_crs(SOURCE_CRS)
    classes: list[str] = []
    component_ids: list[int] = []
    component_node_counts: list[int] = []
    component_line_counts: list[int] = []
    contracted_a: list[str] = []
    contracted_b: list[str] = []
    is_source_component: list[bool] = []

    for _, row in rows.iterrows():
        u = _clean_node(row.NodeA)
        v = _clean_node(row.NodeB)
        cu = SOURCE_NODE if u in sources else u
        cv = SOURCE_NODE if v in sources else v
        contracted_a.append(cu)
        contracted_b.append(cv)

        if cu == cv:
            classes.append("source_self_loop_after_contraction")
            component_ids.append(-1)
            component_node_counts.append(0)
            component_line_counts.append(0)
            is_source_component.append(False)
            continue

        edge = _edge_key(cu, cv)
        if edge in bridges:
            classes.append("bridge")
            component_ids.append(-1)
            component_node_counts.append(0)
            component_line_counts.append(0)
            is_source_component.append(False)
            continue

        component_id = int(node_to_component[cu])
        if component_id != node_to_component[cv]:
            classes.append("between_components_after_bridge_removal")
            component_ids.append(-1)
            component_node_counts.append(0)
            component_line_counts.append(0)
            is_source_component.append(False)
            continue

        classes.append("2_edge_connected_component")
        component_ids.append(component_id)
        component_node_counts.append(component_sizes[component_id])
        component_line_counts.append(component_edge_counts[component_id])
        is_source_component.append(SOURCE_NODE in components[component_id])

    rows["node_a"] = rows["NodeA"].map(_clean_node)
    rows["node_b"] = rows["NodeB"].map(_clean_node)
    rows["contracted_a"] = contracted_a
    rows["contracted_b"] = contracted_b
    rows["connectivity_class"] = classes
    rows["component_id"] = component_ids
    rows["component_nodes"] = component_node_counts
    rows["component_edges"] = component_line_counts
    rows["is_source_component"] = is_source_component
    rows["is_open_tie"] = rows["Status"].astype(str).str.strip() == "0"
    rows["length_m"] = rows.geometry.length.astype(float)

    bridge_gdf = rows[rows["connectivity_class"] == "bridge"].copy()
    component_gdf = rows[rows["connectivity_class"] == "2_edge_connected_component"].copy()
    keep = [
        "node_a",
        "node_b",
        "contracted_a",
        "contracted_b",
        "connectivity_class",
        "component_id",
        "component_nodes",
        "component_edges",
        "is_source_component",
        "is_open_tie",
        "Status",
        "NomV",
        "Feeder",
        "Subest",
        "length_m",
        "geometry",
    ]
    bridge_gdf = bridge_gdf[keep].to_crs("EPSG:3857")
    component_gdf = component_gdf[keep].to_crs("EPSG:3857")
    return bridge_gdf, component_gdf


def write_layers(bridge_gdf: gpd.GeoDataFrame, component_gdf: gpd.GeoDataFrame) -> None:
    pyogrio.write_dataframe(bridge_gdf, GPKG_PATH, layer=BRIDGE_LAYER, driver="GPKG")
    pyogrio.write_dataframe(component_gdf, GPKG_PATH, layer=COMPONENT_LAYER, driver="GPKG")


def update_qgs_project(*, component_ids: list[int]) -> None:
    text = QGS_PATH.read_text(encoding="utf-8")
    for layer_id in (BRIDGE_LAYER_ID, COMPONENT_LAYER_ID):
        text = _remove_existing_layer_references(text, layer_id)

    bridge_tree = _layer_tree_xml(
        layer_id=BRIDGE_LAYER_ID,
        name="MV bridge lines",
        layer=BRIDGE_LAYER,
    )
    component_tree = _layer_tree_xml(
        layer_id=COMPONENT_LAYER_ID,
        name="MV 2-edge-connected component lines",
        layer=COMPONENT_LAYER,
    )
    text = _insert_after_layer_tree_id(text, TIE_LAYER_ID, bridge_tree + "\n" + component_tree)

    bridge_item = f"      <item>{BRIDGE_LAYER_ID}</item>"
    component_item = f"      <item>{COMPONENT_LAYER_ID}</item>"
    text = _insert_after_custom_order_id(text, TIE_LAYER_ID, bridge_item + "\n" + component_item)

    bridge_maplayer = _single_line_maplayer_xml(
        layer_id=BRIDGE_LAYER_ID,
        name="MV bridge lines",
        layer=BRIDGE_LAYER,
        color="55,55,55,255",
        width="0.75",
        style="solid",
    )
    component_maplayer = _categorized_component_maplayer_xml(
        layer_id=COMPONENT_LAYER_ID,
        name="MV 2-edge-connected component lines",
        layer=COMPONENT_LAYER,
        component_ids=component_ids,
    )
    text = text.replace(
        "  </projectlayers>",
        bridge_maplayer + "\n" + component_maplayer + "\n  </projectlayers>",
    )
    QGS_PATH.write_text(text, encoding="utf-8")


def _remove_existing_layer_references(text: str, layer_id: str) -> str:
    text = re.sub(rf"\n\s*<layer-tree-layer[^\n]*id=\"{re.escape(layer_id)}\"[^\n]*/>", "", text)
    text = re.sub(rf"\n\s*<item>{re.escape(layer_id)}</item>", "", text)
    text = re.sub(
        rf"\n\s*<maplayer\b(?:(?!</maplayer>).)*?<id>{re.escape(layer_id)}</id>(?:(?!</maplayer>).)*?</maplayer>",
        "",
        text,
        flags=re.DOTALL,
    )
    return text


def _insert_after_layer_tree_id(text: str, layer_id: str, insertion: str) -> str:
    pattern = rf"(^\s*<layer-tree-layer[^\n]*id=\"{re.escape(layer_id)}\"[^\n]*/>)"
    return re.sub(pattern, rf"\1\n{insertion}", text, count=1, flags=re.MULTILINE)


def _insert_after_custom_order_id(text: str, layer_id: str, insertion: str) -> str:
    pattern = rf"(^\s*<item>{re.escape(layer_id)}</item>)"
    return re.sub(pattern, rf"\1\n{insertion}", text, count=1, flags=re.MULTILINE)


def _layer_tree_xml(*, layer_id: str, name: str, layer: str) -> str:
    return (
        f'    <layer-tree-layer checked="Qt::Unchecked" expanded="1" id="{layer_id}" '
        f'name="{html.escape(name)}" providerKey="ogr" '
        f'source="./P2U_MV_physical_network_3857.gpkg|layername={layer}"/>'
    )


def _single_line_maplayer_xml(
    *,
    layer_id: str,
    name: str,
    layer: str,
    color: str,
    width: str,
    style: str,
) -> str:
    return f"""    <maplayer type="vector" hasScaleBasedVisibilityFlag="0" simplifyDrawingHints="1" simplifyAlgorithm="0" simplifyDrawingTol="1" simplifyLocal="1" labelsEnabled="0" styleCategories="AllStyleCategories">
      <id>{layer_id}</id>
      <datasource>./P2U_MV_physical_network_3857.gpkg|layername={layer}</datasource>
      <keywordList/>
      <layername>{html.escape(name)}</layername>
{_srs_xml()}
      <provider encoding="UTF-8">ogr</provider>
      <renderer-v2 type="singleSymbol" enableorderby="0" forceraster="0" symbollevels="0">
        <symbols>
{_line_symbol_xml("0", color=color, width=width, style=style)}
        </symbols>
      </renderer-v2>
      <geometryOptions removeDuplicateNodes="0" geometryPrecision="0">
        <activeChecks/>
        <checkConfiguration/>
      </geometryOptions>
    </maplayer>"""


def _categorized_component_maplayer_xml(
    *,
    layer_id: str,
    name: str,
    layer: str,
    component_ids: list[int],
) -> str:
    palette = [
        "196,72,72,255",
        "70,132,190,255",
        "166,105,42,255",
        "132,88,175,255",
        "42,145,120,255",
        "210,145,38,255",
        "90,90,90,255",
        "190,85,145,255",
        "92,118,55,255",
        "120,104,205,255",
    ]
    categories = []
    symbols = []
    for symbol_idx, component_id in enumerate(component_ids):
        color = palette[symbol_idx % len(palette)]
        categories.append(
            f'          <category render="true" symbol="{symbol_idx}" '
            f'value="{component_id}" label="component {component_id}"/>'
        )
        symbols.append(_line_symbol_xml(str(symbol_idx), color=color, width="0.62", style="solid"))
    return f"""    <maplayer type="vector" hasScaleBasedVisibilityFlag="0" simplifyDrawingHints="1" simplifyAlgorithm="0" simplifyDrawingTol="1" simplifyLocal="1" labelsEnabled="0" styleCategories="AllStyleCategories">
      <id>{layer_id}</id>
      <datasource>./P2U_MV_physical_network_3857.gpkg|layername={layer}</datasource>
      <keywordList/>
      <layername>{html.escape(name)}</layername>
{_srs_xml()}
      <provider encoding="UTF-8">ogr</provider>
      <renderer-v2 type="categorizedSymbol" attr="component_id" enableorderby="0" forceraster="0" symbollevels="0">
        <categories>
{chr(10).join(categories)}
        </categories>
        <symbols>
{chr(10).join(symbols)}
        </symbols>
      </renderer-v2>
      <geometryOptions removeDuplicateNodes="0" geometryPrecision="0">
        <activeChecks/>
        <checkConfiguration/>
      </geometryOptions>
    </maplayer>"""


def _line_symbol_xml(name: str, *, color: str, width: str, style: str) -> str:
    return f"""          <symbol name="{html.escape(name)}" type="line" alpha="1" clip_to_extent="1" force_rhr="0">
            <layer class="SimpleLine" enabled="1" locked="0" pass="0">
              <Option type="Map">
                <Option name="line_color" type="QString" value="{color}"/>
                <Option name="line_style" type="QString" value="{style}"/>
                <Option name="line_width" type="QString" value="{width}"/>
                <Option name="line_width_unit" type="QString" value="MM"/>
              </Option>
            </layer>
          </symbol>"""


def _srs_xml() -> str:
    return """      <srs>
      <spatialrefsys>
        <wkt></wkt>
        <proj4>+proj=merc +a=6378137 +b=6378137 +lat_ts=0 +lon_0=0 +x_0=0 +y_0=0 +k=1 +units=m +nadgrids=@null +wktext +no_defs</proj4>
        <srsid>3857</srsid>
        <srid>3857</srid>
        <authid>EPSG:3857</authid>
        <description>WGS 84 / Pseudo-Mercator</description>
        <projectionacronym>merc</projectionacronym>
        <ellipsoidacronym>EPSG:7030</ellipsoidacronym>
        <geographicflag>false</geographicflag>
      </spatialrefsys>
    </srs>"""


def update_readme(n_bridges: int, n_component_lines: int, n_components: int) -> None:
    text = README_PATH.read_text(encoding="utf-8")
    addition = f"""\

## Optional Connectivity Layers

- `{BRIDGE_LAYER}`: MV physical lines that are bridges after contracting all HVMV source nodes into one reliability source. Features: `{n_bridges}`.
- `{COMPONENT_LAYER}`: non-bridge MV physical lines labeled by `component_id` for the 2-edge-connected component after source contraction. Features: `{n_component_lines}` across `{n_components}` components.

These layers are unchecked by default in the QGIS project. Turn them on above the base MV line layers to inspect where the physical network is tree-like versus 2-edge-connected.
"""
    if "## Optional Connectivity Layers" in text:
        text = re.sub(r"\n## Optional Connectivity Layers\n.*", addition, text, flags=re.DOTALL)
    else:
        text = text.rstrip() + "\n" + addition
    README_PATH.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
