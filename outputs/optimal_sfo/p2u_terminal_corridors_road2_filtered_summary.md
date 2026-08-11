# P2U Terminal Corridor Preparation

This is a temporary inspection artifact for the street-constrained optimization pipeline.
Street junctions are used only to route corridor geometry. The corridor graph nodes are transformer/source terminals.

## Outputs

- QGIS project: `p2u_terminal_corridors_road2_filtered.qgs`
- GeoPackage: `p2u_terminal_corridors_road2_filtered_3857.gpkg`

## Counts

- Full street graph nodes: `71540`
- Full street graph edges: `81206`
- Source-side road 2-edge nodes kept before terminal graph: `58396`
- Source-side road 2-edge edges kept before terminal graph: `66860`
- Terminal nodes before source-side 2-edge selection: `9529`
- Terminal edges before source-side 2-edge selection: `74322`
- Direct road edges between terminals: `1659`
- Non-terminal road components with at least two terminal boundaries: `5940`
- Max boundary terminals on one non-terminal road component: `196`
- Selected source-side terminal 2-edge nodes: `9529`
- Selected source-side terminal 2-edge edges: `74322`
- Transformers in selected component: `10122` / `12373`
- Sources in selected component: `15` / `15`
- Terminal road nodes: `9529`
- Transformer terminal road nodes: `9518`
- Source terminal road nodes: `15`
- Combined source-transformer terminal road nodes: `4`
- Raw terminal-minor edges in selected component: `74322`
- Transformer-transformer corridor edges: `74065`
- Corridor edges touching a source terminal: `257`
- Corridor edges between transformer-only terminals: `74065`
- Mean raw terminal-minor edge length: `1296.82` m
- Max raw terminal-minor edge length: `5467.79` m
- ILP nodes after degree-2 terminal-chain contraction: `5310`
- ILP transformer nodes: `5299`
- ILP source nodes: `15`
- ILP edges: `70103`
- Mean ILP edge length: `1374.87` m
- Max ILP edge length: `5467.79` m
- Internal transformers moved onto ILP edges: `4507`
- Internal transformer capacity moved onto ILP edges: `334575.0` kVA
- Transformer capacity in selected component: `879720.0` kVA

## Method

1. Attach transformers and HVMV sources to nearest street nodes.
2. Remove terminal nodes from the road graph.
3. For every connected component of non-terminal road nodes, connect all boundary terminals by shortest road-path distance through that component.
4. Add direct terminal-terminal road edges where they exist.
5. Compute the source-side 2-edge-connected component on this terminal-minor graph after contracting source terminals.
6. Contract degree-2 transformer/source-terminal chains into `ilp_edges`, preserving internal transformer count/capacity on the edge.

The resulting graph is intended for visual validation before the ILP stage.
