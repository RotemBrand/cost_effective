# P2U Old-New Reliability Comparison

- Mean edge failure probability target: `0.0005`
- Generalized structural chains computed: `True`

Generalized structural chains were computed with the selected 3-edge component method.

| Metric | Original P2U MV | Optimized backbone + trees | New - old |
|---|---:|---:|---:|
| Source-contracted R | `171` | `50` | `-121` |
| Source-contracted bridges | `8921` | `5206` | `-3715` |
| Raw nodes | `12572` | `10516` | `-2056` |
| Raw edges | `12728` | `10551` | `-2177` |
| Total length km | `589.339` | `587.791` | `-1.54817` |
| Reliability weight kW | `661467` | `661362` | `-105.22` |
| 2-edge components | `8922` | `5207` | `-3715` |
| Structure graph nodes | `1426` | `1736` | `310` |
| Structure graph edges | `1590` | `1782` | `192` |
| Regular chains | `1590` | `1785` | `195` |
| Float total risk | `0.0022664` | `0.000674192` | `-0.00159221` |
| Float tree risk | `0.00226025` | `0.000341193` | `-0.00191906` |
| Float nonbridge section risk | `0` | `6.76344e-05` | `6.76344e-05` |
| Float internal regular-chain risk | `3.6127e-07` | `6.38796e-06` | `6.02669e-06` |
| Float structural generalized-chain risk | `5.79277e-06` | `0.000258977` | `0.000253185` |
| Poly total p1 coeff | `4.54049` | `0.817909` | `-3.72258` |
| Poly total p2 coeff | `-15.5583` | `1065.27` | `1080.82` |
| Poly value at p_mean | `0.00226635` | `0.000675271` | `-0.00159108` |

## Decomposition Meaning

- `tree`: bridge / 1-connected contribution after source contraction.
- `nonbridge_section`: first-order risk from load stored on non-bridge section edges.
- `internal_regular_chains`: second-order risk from two cuts inside regular chains.
- `structural_generalized_chains`: second-order generalized-chain structural term. It is zero when `--generalized` is not used.

## Runtime

- Original decomposition: `6.22` s
- New decomposition: `8.09` s
