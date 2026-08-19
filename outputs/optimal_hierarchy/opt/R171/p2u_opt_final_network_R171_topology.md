# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `11842` / `11998`
- Source-contracted nodes/edges: `11828` / `11998`
- Source-contracted cycle rank R: `171`
- Source-contracted bridge edges: `6532`
- Backbone edges: `5466`
- Tree attachment edges: `6532`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `12373`
- Contracted edge-load transformers: `0`
- Transformers on backbone, including edge load: `5615`
- Tree transformers: `6758`
- Percent transformers on backbone: `45.381%`
- Capacity on backbone: `545145.0` kVA / `1100120.0` kVA
- Demand on backbone: `379153.260` kW / `661466.740` kW
- Percent demand on backbone: `57.320%`

## Length

- Backbone length: `691976.508` m
- Tree physical road-union length: `216849.203` m
- Total physical length: `908825.712` m
- Terminal attachment distance sum: `216849.203` m
- Mean/max terminal distance to backbone: `33.198` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `6.6008256e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
