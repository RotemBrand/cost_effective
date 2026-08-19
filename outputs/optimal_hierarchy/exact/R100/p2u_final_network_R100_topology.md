# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `10543` / `10628`
- Source-contracted nodes/edges: `10529` / `10628`
- Source-contracted cycle rank R: `100`
- Source-contracted bridge edges: `5233`
- Backbone edges: `5395`
- Tree attachment edges: `5233`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `10989`
- Contracted edge-load transformers: `1384`
- Transformers on backbone, including edge load: `6999`
- Tree transformers: `5374`
- Percent transformers on backbone: `56.567%`
- Capacity on backbone: `650620.0` kVA / `1100120.0` kVA
- Demand on backbone: `427511.510` kW / `661466.740` kW
- Percent demand on backbone: `64.631%`

## Length

- Backbone length: `399755.967` m
- Tree physical road-union length: `188372.003` m
- Total physical length: `588127.970` m
- Terminal attachment distance sum: `188372.003` m
- Mean/max terminal distance to backbone: `35.997` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `9.0354485e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
