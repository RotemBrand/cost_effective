# P2U Final Optimized Network Analysis

The analyzed graph is a switch-section graph. Transformer chains are represented as edge load; street-forest versions may also include zero-load branch nodes.

## Topology

- Final graph nodes/edges: `11841` / `11890`
- Source-contracted nodes/edges: `11831` / `11890`
- Source-contracted cycle rank R: `60`
- Source-contracted bridge edges: `6531`
- Backbone edges: `5359`
- Tree attachment edges: `6531`
- Sources: `11`

## Transformers And Load

- Represented transformers: `12369` / `12373`
- Node transformers: `12369`
- Contracted edge-load transformers: `0`
- Transformers on backbone, including edge load: `5611`
- Tree transformers: `6758`
- Percent transformers on backbone: `45.363%`
- Capacity on backbone: `544870.0` kVA / `1099845.0` kVA
- Demand on backbone: `379048.040` kW / `661361.520` kW
- Percent demand on backbone: `57.313%`

## Length

- Backbone length: `562812.236` m
- Tree physical road-union length: `216841.423` m
- Total physical length: `779653.659` m
- Terminal attachment distance sum: `216841.423` m
- Mean/max terminal distance to backbone: `33.202` m / `816.951` m

## Failure Probabilities

- Target mean edge probability: `0.0005`
- Fitted failure rate per meter: `7.6251807e-06`

## Reliability

- Monte Carlo: `not run`

## Notes

- The tree attachment edge length is the road shortest-path distance to the selected backbone.
- The physical tree length uses the union of road edges, so it avoids double-counting shared attachment paths.
