# SFO Random-Switch Original-Style Figure Stats

Topology source: saved paper networks from `data/real_networks.nxjson`.

Switch layer: non-tree edges are marked as tie switches; remaining switches are placed randomly on closed edges until the total switch count matches the corresponding SMART-DS feeder.

MCMC settings: `p=0.0005` mean section-edge failure probability for each original switched graph, `T_days=1825`, `mean_cycle_days=0.01`, simulation seed `100`, switch seed `2026`.

| name | N | M | R | R_optimized | reference_switches | original_switches | optimized_switches | original_ties | optimized_ties | total_weight | optimal_network_weight | saidi | optimal_network_saidi | weight_ratio | R_ratio | rel_ratio | edge_failure_rate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SFO Pacific | 944 | 951 | 8 | 8 | 737 | 737 | 737 | 8 | 8 | 56049.4 | 61227.7 | 0.00159261 | 0.000375989 | -0.0923869 | 0 | 0.763916 | 8.1892e-06 |
| SFO Davidson | 489 | 494 | 6 | 6 | 377 | 377 | 377 | 6 | 6 | 40381.9 | 43955.3 | 0.00130613 | 0.000195793 | -0.0884893 | 0 | 0.850096 | 5.95564e-06 |

Index definitions:

- `Z_w = weight_ratio = 1 - W_optimized / W_original`.
- `Z_R = R_ratio = 1 - R_optimized / R_original`.
- `Z_F = rel_ratio = 1 - F_optimized / F_original`.
- `F` is computed after random switch placement and switch-section contraction.
