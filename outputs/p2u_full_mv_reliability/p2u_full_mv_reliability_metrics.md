# P2U Full MV Reliability Metrics

## Topology

- Nodes: `12572`
- Unique MV line edges: `12728`
- Raw MV line features: `12732`
- Normally open tie edges R: `171`
- Physical cycle-rank total ties: `157`
- MV switching-device features: `16284`
- Closed operating components: `15`
- Physical components with ties: `1`

## 2-Edge-Connected Backbone

- Source 2-edge-connected backbone nodes: `3562`
- Nodes in any nontrivial 2-edge-connected component: `3642`
- Largest 2-edge-connected component nodes: `3562`
- Bridge edges after source contraction: `8921`

## Consumers On Source 2-Edge-Connected Backbone

- Transformer load nodes: `12373`
- Transformer load nodes 2-edge-connected to source: `3509`
- Percent transformer load nodes 2-edge-connected to source: `28.360%`
- Demand 2-edge-connected to source: `234948.610` kW / `661466.740` kW
- Percent demand 2-edge-connected to source: `35.519%`
- Customer count 2-edge-connected to source: `60317.000` / `191140.000`
- Percent customers 2-edge-connected to source: `31.556%`

## LV Assignment

- LV consumer points assigned to MV transformers: `46379` / `46379`
- Consumer points assigned by nearest-transformer fallback: `196`
- Nearest-transformer fallback mean/max distance: `44.60` m / `203.92` m
- Assigned demand: `661466.740` kW / `661466.740` kW
- Transformer nodes with assigned consumers: `12373`

## SAIDI Monte Carlo

- Method: `stationary_independent_edge_state_monte_carlo`
- Edge failure probability p: `0.0005`
- Samples: `50000`
- Seed: `20260810`
- Mean disconnected load fraction SAIDI: `0.0019909671`
- 95% CI half-width: `2.45e-05`
- Expected failed edges per sample: `6.337`
- Mean failed edges per sample: `6.337`
- Runtime seconds: `775.10`

## Assumptions

- MV network uses all `NomV` between 1 kV and 40 kV.
- Normally open ties are `Status = 0` MV lines.
- SAIDI is a stationary independent-edge Monte Carlo estimate, not the slower event-loop simulator.
- LV consumer demand is aggregated to MV transformer nodes by LV connected component.
- LV consumers not found in an LV line component are assigned to the nearest distribution transformer.
- 2-connected means 2-edge-connected relative to the contracted source set, because the failure model is edge/line failure.