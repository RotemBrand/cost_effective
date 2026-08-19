# P2U OPT Hierarchy: try_001_opt_equal_chain

## Results

This experiment replaces the road-corridor ILP backbone with the OPT equal-chain construction on the same backbone-terminal set. The final graph still uses the street-forest attachment for non-backbone transformers.

The OPT backbone is Euclidean, so the backbone cost is not road-constrained. The tree attachment cost remains street-based.

| network | r_request | length_km | r_theory | z_w | z_r | z_f | z_f_p | z_f_p2 | risk_total | risk_o_p | risk_o_p2 | risk_tree | risk_section | risk_internal | risk_structural | bridges | generalized_chains | gen_lambda_mean_km | gen_lambda_sigma_over_mean | gen_lambda_max_km | ilp_runtime_s | decomposition_runtime_s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Original P2U MV |  | 5.89e2 | 171 | 0 | 0 | 0 | 0 | 0 | 2.27e-3 | 2.26e-3 | 6.15e-6 | 2.26e-3 | 0 | 3.61e-7 | 5.79e-6 | 8921 | 324 | 3.24e-1 | 1.45e0 | 4.04e0 |  | 7.72e0 |
| OPT hierarchy R50 | 50 | 7.80e2 | 60 | -3.23e-1 | 6.49e-1 | 8.28e-1 | 8.64e-1 | -1.22e1 | 3.89e-4 | 3.08e-4 | 8.11e-5 | 3.08e-4 | 0 | 3.97e-6 | 7.71e-5 | 6531 | 158 | 2.37e0 | 7.10e-1 | 1.04e1 | 5.60e1 | 2.75e1 |

## Algorithm

1. Load the same terminal set used by the road-corridor ILP backbone.
2. Run balanced chain clustering and OPT structure construction for the requested `R`.
3. Export the OPT graph as a backbone solution with straight-line Euclidean edges.
4. Reuse the final-network builder to attach all remaining transformers as a street forest.
5. Run the deterministic switch-aware reliability decomposition with generalized chains.

## Implementation

- Script: `figures/optimal_hierarchy/opt_algorithm_sweep.py`
- Output table: `outputs\optimal_hierarchy\opt\try_001_opt_equal_chain\p2u_opt_hierarchy_table.csv`
- Output JSON: `outputs\optimal_hierarchy\opt\try_001_opt_equal_chain\p2u_opt_hierarchy_table.json`

Parameters:

```json
{
  "p_mean": 0.0005,
  "seed": 20260812,
  "kmeans_max_iter": 3,
  "strc_n_init_iters": 1,
  "strc_exact_vertices": true,
  "strc_trip_nearest_vertices": 24,
  "chain_n_init_iters": 1,
  "local_fix_max_changes": 0,
  "local_fix_max_risk_gain": 0.2,
  "generalized_method": "projection",
  "tree_mode": "street_forest"
}
```

Stage log:

```json
{
  "parameters": {
    "p_mean": 0.0005,
    "seed": 20260812,
    "kmeans_max_iter": 3,
    "strc_n_init_iters": 1,
    "strc_exact_vertices": true,
    "strc_trip_nearest_vertices": 24,
    "chain_n_init_iters": 1,
    "local_fix_max_changes": 0,
    "local_fix_max_risk_gain": 0.2,
    "generalized_method": "projection",
    "tree_mode": "street_forest"
  },
  "corridor_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_terminal_corridors_road2_k10_3857.gpkg",
  "backbone_terminal_count": 5310,
  "cases": {
    "50": {
      "backbone": {
        "algorithm": "OPT_equal_chain_backbone",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\opt\\try_001_opt_equal_chain\\R050\\p2u_opt_backbone_R50_3857.gpkg",
        "redundancy_constraint": 50,
        "max_redundancy_constraint": null,
        "status_name": "OPT_CONSTRUCTED",
        "stop_reason": "opt_algorithm_completed",
        "runtime_s": 56.02468799997587,
        "seed": 20260812,
        "kmeans_max_iter": 3,
        "strc_n_init_iters": 1,
        "strc_exact_vertices": true,
        "strc_trip_nearest_vertices": 24,
        "chain_n_init_iters": 1,
        "local_fix": {
          "enabled": false
        },
        "input_nodes": 5310,
        "solution_nodes": 5310,
        "solution_edges": 5359,
        "solution_cycle_rank": 50,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "solution_transformer_nodes": 5295,
        "solution_source_nodes": 11,
        "objective_length_m": 562812.2363966871,
        "mean_edge_length_m": 105.02187654351317,
        "max_edge_length_m": 2131.4006322548075,
        "tie_edges": 50,
        "tree_mode": "street_forest"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\opt\\try_001_opt_equal_chain\\R050\\p2u_opt_final_network_R50_risk.json",
      "case_runtime_s": 42.477153199957684
    }
  },
  "total_runtime_s": 42.655042900005355
}
```
