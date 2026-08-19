# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `10516` / `10591`
- Source-contracted nodes/edges: `10502` / `10591`
- Source-contracted cycle rank R: `90`
- Source-contracted bridge edges: `5206`
- Backbone edges: `5385`
- Tree attachment edges: `5206`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `10960`
- Contracted edge-load transformers: `1413`
- Transformers on backbone, including edge load: `7028`
- Tree transformers: `5345`
- Percent transformers on backbone: `56.801%`
- Capacity on backbone: `651470.0` kVA / `1100120.0` kVA
- Demand on backbone: `427891.290` kW / `661466.740` kW
- Percent demand on backbone: `64.688%`

## Length

- Backbone length: `400812.272` m
- Tree physical road-union length: `187749.290` m
- Total physical length: `588561.562` m
- Terminal attachment distance sum: `187749.290` m
- Mean/max terminal distance to backbone: `36.064` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `8.9973596e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
