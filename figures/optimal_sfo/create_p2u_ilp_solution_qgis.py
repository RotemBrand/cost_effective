from __future__ import annotations

import html
import json
from pathlib import Path

import geopandas as gpd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo"
SUMMARY_JSON = OUTPUT_DIR / "p2u_ilp_2edge_solution_summary.json"
QGS_PATH = OUTPUT_DIR / "p2u_ilp_2edge_solution.qgs"


def _line_symbol(color: str, width: str) -> str:
    return f"""      <renderer-v2 type="singleSymbol">
        <symbols>
          <symbol alpha="1" type="line" name="0">
            <layer enabled="1" class="SimpleLine">
              <Option type="Map">
                <Option name="line_color" type="QString" value="{color}"/>
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


def main() -> None:
    summary = json.loads(SUMMARY_JSON.read_text(encoding="utf-8"))
    gpkg = Path(summary["output_gpkg"])
    gpkg_name = gpkg.name

    # Read once to fail loudly if the expected layers are missing.
    gpd.read_file(gpkg, layer="solution_edges")
    gpd.read_file(gpkg, layer="solution_transformer_nodes")
    gpd.read_file(gpkg, layer="solution_source_nodes")

    layers = [
        ("solution_edges", "ILP solution edges", "solution_edges_eb3db9ef", "Line", _line_symbol("174,56,132,255", "0.65")),
        ("solution_transformer_nodes", "ILP transformer nodes", "solution_transformer_nodes_8e663113", "Point", _point_symbol("255,255,255,0", "203,64,60,255", "1.2", False)),
        ("solution_source_nodes", "ILP source nodes", "solution_source_nodes_1bd69f03", "Point", _point_symbol("18,18,18,255", "255,255,255,255", "2.6", True)),
    ]
    tree = "\n".join(
        f'      <layer-tree-layer checked="1" id="{layer_id}" name="{html.escape(name)}" source="./{gpkg_name}|layername={layer}"/>'
        for layer, name, layer_id, _, _ in layers
    )
    project_layers = "\n".join(
        _map_layer(gpkg_name, layer, name, layer_id, geometry, renderer)
        for layer, name, layer_id, geometry, renderer in layers
    )
    qgs = f"""<!DOCTYPE qgis PUBLIC 'http://mrcc.com/qgis.dtd' 'SYSTEM'>
<qgis version="3.34.0" projectname="P2U ILP 2-edge solution">
  <homePath path="."/>
  <title>P2U ILP 2-edge solution</title>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">
    <customproperties/>
    <layer-tree-group checked="Qt::Checked" expanded="1" name="P2U ILP 2-edge solution">
{tree}
    </layer-tree-group>
  </layer-tree-group>
  <projectlayers>
{project_layers}
  </projectlayers>
</qgis>
"""
    QGS_PATH.write_text(qgs, encoding="utf-8")
    print(QGS_PATH.name)


if __name__ == "__main__":
    main()
