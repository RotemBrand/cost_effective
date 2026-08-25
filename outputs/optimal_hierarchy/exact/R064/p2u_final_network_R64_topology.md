# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `9803` / `9852`
- Source-contracted nodes/edges: `9789` / `9852`
- Source-contracted cycle rank R: `64`
- Source-contracted bridge edges: `4493`
- Backbone edges: `5359`
- Tree attachment edges: `4493`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `10223`
- Contracted edge-load transformers: `2150`
- Transformers on backbone, including edge load: `7765`
- Tree transformers: `4608`
- Percent transformers on backbone: `62.758%`
- Capacity on backbone: `710845.0` kVA / `1100120.0` kVA
- Demand on backbone: `453906.590` kW / `661466.740` kW
- Percent demand on backbone: `68.621%`

## Length

- Backbone length: `442855.259` m
- Tree physical road-union length: `171994.857` m
- Total physical length: `614850.116` m
- Terminal attachment distance sum: `171994.857` m
- Mean/max terminal distance to backbone: `38.281` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `8.0117087e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
