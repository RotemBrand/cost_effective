# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `10534` / `10669`
- Source-contracted nodes/edges: `10520` / `10669`
- Source-contracted cycle rank R: `150`
- Source-contracted bridge edges: `5224`
- Backbone edges: `5445`
- Tree attachment edges: `5224`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `10974`
- Contracted edge-load transformers: `1399`
- Transformers on backbone, including edge load: `7014`
- Tree transformers: `5359`
- Percent transformers on backbone: `56.688%`
- Capacity on backbone: `651620.0` kVA / `1100120.0` kVA
- Demand on backbone: `427859.230` kW / `661466.740` kW
- Percent demand on backbone: `64.683%`

## Length

- Backbone length: `399974.446` m
- Tree physical road-union length: `188359.007` m
- Total physical length: `588333.454` m
- Terminal attachment distance sum: `188359.007` m
- Mean/max terminal distance to backbone: `36.056` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `9.067137e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
