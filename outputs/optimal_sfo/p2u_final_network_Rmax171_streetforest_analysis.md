# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `10504` / `10656`
- Source-contracted nodes/edges: `10490` / `10656`
- Source-contracted cycle rank R: `167`
- Source-contracted bridge edges: `5194`
- Backbone edges: `5462`
- Tree attachment edges: `5194`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `10943`
- Contracted edge-load transformers: `1430`
- Transformers on backbone, including edge load: `7045`
- Tree transformers: `5328`
- Percent transformers on backbone: `56.938%`
- Capacity on backbone: `654120.0` kVA / `1100120.0` kVA
- Demand on backbone: `428917.810` kW / `661466.740` kW
- Percent demand on backbone: `64.843%`

## Length

- Backbone length: `403139.909` m
- Tree physical road-union length: `187726.224` m
- Total physical length: `590866.133` m
- Terminal attachment distance sum: `187726.224` m
- Mean/max terminal distance to backbone: `36.143` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `9.0172709e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
