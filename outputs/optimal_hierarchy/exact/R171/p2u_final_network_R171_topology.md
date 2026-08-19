# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `10531` / `10687`
- Source-contracted nodes/edges: `10517` / `10687`
- Source-contracted cycle rank R: `171`
- Source-contracted bridge edges: `5221`
- Backbone edges: `5466`
- Tree attachment edges: `5221`
- Sources: `15`

## Transformers And Load

- Represented transformers: `12373` / `12373`
- Node transformers: `10971`
- Contracted edge-load transformers: `1402`
- Transformers on backbone, including edge load: `7017`
- Tree transformers: `5356`
- Percent transformers on backbone: `56.712%`
- Capacity on backbone: `652020.0` kVA / `1100120.0` kVA
- Demand on backbone: `428006.200` kW / `661466.740` kW
- Percent demand on backbone: `64.706%`

## Length

- Backbone length: `400349.125` m
- Tree physical road-union length: `188303.665` m
- Total physical length: `588652.790` m
- Terminal attachment distance sum: `188303.665` m
- Mean/max terminal distance to backbone: `36.067` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `9.0775073e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
