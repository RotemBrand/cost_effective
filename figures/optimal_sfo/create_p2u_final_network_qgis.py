from __future__ import annotations

import argparse
import html
import json
import sys
from pathlib import Path

import geopandas as gpd
import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.p2u_final_network import (
    FINAL_NETWORK_GPKG,
    FINAL_NETWORK_METADATA,
    OUTPUT_DIR,
)


DEFAULT_QGS = OUTPUT_DIR / "p2u_final_network.qgs"


def _line_symbol(color: str, width: str, style: str = "solid", capstyle: str = "round") -> str:
    return f"""      <renderer-v2 type="singleSymbol">
        <symbols>
          <symbol alpha="1" type="line" name="0">
            <layer enabled="1" class="SimpleLine">
              <Option type="Map">
                <Option name="line_color" type="QString" value="{color}"/>
                <Option name="line_style" type="QString" value="{style}"/>
                <Option name="line_width" type="QString" value="{width}"/>
                <Option name="line_width_unit" type="QString" value="MM"/>
                <Option name="capstyle" type="QString" value="{capstyle}"/>
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
                <Option name="outline_width" type="QString" value="0.35"/>
                <Option name="outline_width_unit" type="QString" value="MM"/>
                <Option name="size" type="QString" value="{size}"/>
                <Option name="size_unit" type="QString" value="MM"/>
              </Option>
            </layer>
          </symbol>
        </symbols>
      </renderer-v2>"""


def _map_layer(
    *,
    gpkg_name: str,
    layer: str,
    name: str,
    layer_id: str,
    geometry: str,
    renderer: str,
    checked: str,
) -> tuple[str, str]:
    tree = f'      <layer-tree-layer checked="{checked}" id="{layer_id}" name="{html.escape(name)}" source="./{gpkg_name}|layername={layer}"/>'
    project = f"""  <maplayer type="vector" geometry="{geometry}">
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
    return tree, project


def _write_switch_style_layers(gpkg_path: Path, style_gpkg: Path) -> dict:
    backbone = gpd.read_file(gpkg_path, layer="final_backbone_edges").reset_index(drop=True)
    if backbone.empty:
        raise ValueError(f"{gpkg_path} has no final_backbone_edges rows")

    graph = nx.Graph()
    for idx, row in backbone.iterrows():
        u = str(row.terminal_a)
        v = str(row.terminal_b)
        if u == v:
            continue
        length = float(row.length_m)
        graph.add_edge(u, v, weight=length, row_id=int(idx))

    normally_closed_rows: set[int] = set()
    for component in nx.connected_components(graph):
        subgraph = graph.subgraph(component)
        for _, _, data in nx.minimum_spanning_edges(subgraph, data=True):
            normally_closed_rows.add(int(data["row_id"]))

    normally_open = backbone.loc[~backbone.index.isin(normally_closed_rows)].copy()
    normally_closed = backbone.loc[backbone.index.isin(normally_closed_rows)].copy()
    normally_closed["switch_state"] = "normally_closed"
    normally_open["switch_state"] = "normally_open_tie"

    if style_gpkg.exists():
        style_gpkg.unlink()
    normally_closed.to_file(style_gpkg, layer="backbone_normally_closed", driver="GPKG")
    normally_open.to_file(style_gpkg, layer="backbone_normally_open_ties", driver="GPKG")
    return {
        "style_gpkg": str(style_gpkg),
        "normally_closed_backbone_edges": int(len(normally_closed)),
        "normally_open_tie_edges": int(len(normally_open)),
    }


def write_qgis_project(
    *,
    gpkg_path: Path = FINAL_NETWORK_GPKG,
    qgs_path: Path = DEFAULT_QGS,
) -> dict:
    gpkg_path = Path(gpkg_path)
    qgs_path = Path(qgs_path)
    style_gpkg = qgs_path.with_suffix(".style_layers.gpkg")
    switch_layers = _write_switch_style_layers(gpkg_path, style_gpkg)
    specs = [
        {
            "gpkg_name": style_gpkg.name,
            "layer": "backbone_normally_closed",
            "name": "Backbone normally closed lines",
            "layer_id": "backbone_normally_closed_c142ac02",
            "geometry": "Line",
            "renderer": _line_symbol("203,64,60,255", "0.70"),
            "checked": "1",
        },
        {
            "gpkg_name": style_gpkg.name,
            "layer": "backbone_normally_open_ties",
            "name": "Backbone normally open ties",
            "layer_id": "backbone_normally_open_ties_ad74edf3",
            "geometry": "Line",
            "renderer": _line_symbol("245,157,35,255", "0.95", "dash"),
            "checked": "1",
        },
        {
            "gpkg_name": gpkg_path.name,
            "layer": "final_tree_attachment_edges",
            "name": "Radial tree attachments",
            "layer_id": "final_tree_attachment_edges_a51991ef",
            "geometry": "Line",
            "renderer": _line_symbol("77,175,74,210", "0.35", "dash"),
            "checked": "1",
        },
        {
            "gpkg_name": gpkg_path.name,
            "layer": "final_transformer_nodes",
            "name": "Transformer terminals",
            "layer_id": "final_transformer_nodes_9cc0a42d",
            "geometry": "Point",
            "renderer": _point_symbol("255,255,255,0", "55,55,55,255", "0.95", False),
            "checked": "1",
        },
        {
            "gpkg_name": gpkg_path.name,
            "layer": "final_source_nodes",
            "name": "Sources",
            "layer_id": "final_source_nodes_2d3d5024",
            "geometry": "Point",
            "renderer": _point_symbol("18,18,18,255", "255,255,255,255", "2.5", True),
            "checked": "1",
        },
        {
            "gpkg_name": gpkg_path.name,
            "layer": "final_street_branch_nodes",
            "name": "Street branch nodes",
            "layer_id": "final_street_branch_nodes_7c1db8a2",
            "geometry": "Point",
            "renderer": _point_symbol("255,255,255,0", "120,120,120,180", "0.65", False),
            "checked": "1",
        },
    ]
    trees = []
    projects = []
    for spec in specs:
        tree, project = _map_layer(**spec)
        trees.append(tree)
        projects.append(project)

    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="P2U final optimized network">
  <homePath path="."/>
  <title>P2U final optimized network</title>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="P2U final optimized network">
{chr(10).join(trees)}
    </layer-tree-group>
  </layer-tree-group>
  <projectlayers>
{chr(10).join(projects)}
  </projectlayers>
</qgis>
"""
    qgs_path.write_text(qgs, encoding="utf-8")
    return {"output_qgs": str(qgs_path), "input_gpkg": str(gpkg_path), **switch_layers}


def _gpkg_from_metadata(metadata_json: Path) -> Path:
    if not metadata_json.exists():
        return FINAL_NETWORK_GPKG
    metadata = json.loads(metadata_json.read_text(encoding="utf-8"))
    return Path(metadata.get("output_gpkg", FINAL_NETWORK_GPKG))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a temporary QGIS project for the final P2U optimized network.")
    parser.add_argument("--metadata-json", type=Path, default=FINAL_NETWORK_METADATA)
    parser.add_argument("--gpkg", type=Path, default=None)
    parser.add_argument("--output-qgs", type=Path, default=DEFAULT_QGS)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gpkg = args.gpkg if args.gpkg is not None else _gpkg_from_metadata(args.metadata_json)
    print(json.dumps(write_qgis_project(gpkg_path=gpkg, qgs_path=args.output_qgs), indent=2))


if __name__ == "__main__":
    main()
