# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `11744` / `11829`
- Source-contracted nodes/edges: `11730` / `11829`
- Source-contracted cycle rank R: `100`
- Source-contracted bridge edges: `2215`
- Backbone edges: `9614`
- Tree attachment edges: `2215`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `12373`
- Contracted edge-load transformers: `0`
- Transformers on backbone, including edge load: `10122`
- Tree transformers: `2251`
- Percent transformers on backbone: `81.807%`
- Capacity on backbone: `879720.0` kVA / `1100120.0` kVA
- Demand on backbone: `524630.530` kW / `661466.740` kW
- Percent demand on backbone: `79.313%`

## Length

- Backbone length: `575455.883` m
- Tree physical road-union length: `102055.759` m
- Total physical length: `677511.641` m
- Terminal attachment distance sum: `102055.759` m
- Mean/max terminal distance to backbone: `46.075` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `8.7297393e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
