# P2U Old-New Reliability Comparison

- Mean edge failure probability target: `0.0005`
- Generalized structural chains computed: `True`

Generalized structural chains were computed with the selected 3-edge component method.

| Metric | Original P2U MV | Optimized backbone + trees | New - old |
|---|---:|---:|---:|
| Source-contracted R | `171` | `167` | `-4` |
| Source-contracted bridges | `8921` | `5194` | `-3727` |
| Raw nodes | `12572` | `10504` | `-2068` |
| Raw edges | `12728` | `10656` | `-2072` |
| Total length km | `589.339` | `590.866` | `1.52676` |
| Reliability weight kW | `661467` | `661362` | `-105.22` |
| 2-edge components | `8922` | `5195` | `-3727` |
| Structure graph nodes | `1426` | `1855` | `429` |
| Structure graph edges | `1590` | `2013` | `423` |
| Regular chains | `1590` | `2021` | `431` |
| Float total risk | `0.0022664` | `0.000447646` | `-0.00181876` |
| Float tree risk | `0.00226025` | `0.000342808` | `-0.00191744` |
| Float nonbridge section risk | `0` | `6.73449e-05` | `6.73449e-05` |
| Float internal regular-chain risk | `3.6127e-07` | `4.45952e-06` | `4.09825e-06` |
| Float structural generalized-chain risk | `5.79277e-06` | `3.3034e-05` | `2.72412e-05` |
| Poly total p1 coeff | `4.54049` | `0.820563` | `-3.71993` |
| Poly total p2 coeff | `-15.5583` | `150.021` | `165.579` |
| Poly value at p_mean | `0.00226635` | `0.000447787` | `-0.00181857` |

## Decomposition Meaning

- `tree`: bridge / 1-connected contribution after source contraction.
- `nonbridge_section`: first-order risk from load stored on non-bridge section edges.
- `internal_regular_chains`: second-order risk from two cuts inside regular chains.
- `structural_generalized_chains`: second-order generalized-chain structural term. It is zero when `--generalized` is not used.

## Runtime

- Original decomposition: `7.43` s
- New decomposition: `7.39` s
