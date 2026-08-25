# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `8931` / `8966`
- Source-contracted nodes/edges: `8917` / `8966`
- Source-contracted cycle rank R: `50`
- Source-contracted bridge edges: `3621`
- Backbone edges: `5345`
- Tree attachment edges: `3621`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `9300`
- Contracted edge-load transformers: `3073`
- Transformers on backbone, including edge load: `8688`
- Tree transformers: `3685`
- Percent transformers on backbone: `70.217%`
- Capacity on backbone: `779970.0` kVA / `1100120.0` kVA
- Demand on backbone: `483023.370` kW / `661466.740` kW
- Percent demand on backbone: `73.023%`

## Length

- Backbone length: `514252.832` m
- Tree physical road-union length: `153028.449` m
- Total physical length: `667281.281` m
- Terminal attachment distance sum: `153028.449` m
- Mean/max terminal distance to backbone: `42.261` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `6.7183062e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
