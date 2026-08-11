# P2U Terminal Corridor Preparation

This is a temporary inspection artifact for the street-constrained optimization pipeline.
Street junctions are used only to route corridor geometry. The corridor graph nodes are transformer/source terminals.

## Outputs

- QGIS project: `p2u_terminal_corridors.qgs`
- GeoPackage: `p2u_terminal_corridors_3857.gpkg`

## Counts

- Full street graph nodes: `71540`
- Full street graph edges: `81206`
- Terminal nodes before source-side 2-edge selection: `11555`
- Terminal edges before source-side 2-edge selection: `93192`
- Direct road edges between terminals: `2184`
- Non-terminal road components with at least two terminal boundaries: `6618`
- Max boundary terminals on one non-terminal road component: `226`
- Selected source-side terminal 2-edge nodes: `10429`
- Selected source-side terminal 2-edge edges: `91903`
- Transformers in selected component: `11107` / `12373`
- Sources in selected component: `15` / `15`
- Terminal road nodes: `11555`
- Transformer terminal road nodes: `11544`
- Source terminal road nodes: `15`
- Combined source-transformer terminal road nodes: `4`
- Raw terminal-minor edges in selected component: `91903`
- Transformer-transformer corridor edges: `91611`
- Corridor edges touching a source terminal: `292`
- Corridor edges between transformer-only terminals: `91611`
- Mean raw terminal-minor edge length: `1340.72` m
- Max raw terminal-minor edge length: `5627.39` m
- ILP nodes after degree-2 terminal-chain contraction: `6333`
- ILP transformer nodes: `6322`
- ILP source nodes: `15`
- ILP edges: `87807`
- Mean ILP edge length: `1403.26` m
- Max ILP edge length: `5627.39` m
- Internal transformers moved onto ILP edges: `4381`
- Internal transformer capacity moved onto ILP edges: `324425.0` kVA
- Transformer capacity in selected component: `989970.0` kVA

## Method

1. Attach transformers and HVMV sources to nearest street nodes.
2. Remove terminal nodes from the road graph.
3. For every connected component of non-terminal road nodes, connect all boundary terminals by shortest road-path distance through that component.
4. Add direct terminal-terminal road edges where they exist.
5. Compute the source-side 2-edge-connected component on this terminal-minor graph after contracting source terminals.
6. Contract degree-2 transformer/source-terminal chains into `ilp_edges`, preserving internal transformer count/capacity on the edge.

The resulting graph is intended for visual validation before the ILP stage.
