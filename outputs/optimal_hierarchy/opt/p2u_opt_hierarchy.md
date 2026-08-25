# P2U OPT Hierarchy: opt

## Results

This experiment replaces the road-corridor ILP backbone with the OPT equal-chain construction. By default it uses the strict road-2edge terminal set: every transformer attached to a source-side 2-edge-connected street node is aggregated by road node and included as a mandatory OPT terminal. The final graph still uses the street-forest attachment for road-1edge transformers.

The OPT backbone is Euclidean, so the backbone cost is not road-constrained. The tree attachment cost remains street-based.

When `target_source_contracted_r=true`, the OPT physical redundancy is reduced by `n_sources - 1` so the final source-contracted theory graph has the requested `R`.

| network | r_request | length_km | r_theory | z_w | z_r | z_f | z_f_p | z_f_p2 | risk_total | risk_o_p | risk_o_p2 | risk_tree | risk_section | risk_internal | risk_structural | bridges | generalized_chains | gen_lambda_mean_km | gen_lambda_sigma_over_mean | gen_lambda_max_km | ilp_runtime_s | decomposition_runtime_s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Original P2U MV |  | 5.89e2 | 171 | 0 | 0 | 0 | 0 | 0 | 2.27e-3 | 2.26e-3 | 6.15e-6 | 2.26e-3 | 0 | 3.61e-7 | 5.79e-6 | 8921 | 324 | 3.24e-1 | 1.45e0 | 4.04e0 |  | 8.78e0 |
| OPT hierarchy R50 | 50 | 6.24e2 | 50 | -5.96e-2 | 7.08e-1 | 7.49e-1 | 8.82e-1 | -4.81e1 | 5.70e-4 | 2.68e-4 | 3.02e-4 | 2.68e-4 | 0 | 3.57e-5 | 2.66e-4 | 2215 | 112 | 3.61e0 | 5.96e-1 | 1.20e1 | 1.11e2 | 7.89e0 |
| OPT hierarchy R100 | 100 | 6.78e2 | 100 | -1.50e-1 | 4.15e-1 | 8.46e-1 | 8.82e-1 | -1.23e1 | 3.50e-4 | 2.68e-4 | 8.20e-5 | 2.68e-4 | 0 | 2.76e-5 | 5.44e-5 | 2215 | 214 | 1.71e0 | 7.85e-1 | 1.19e1 | 2.44e2 | 8.35e0 |

## Algorithm

1. Load and aggregate the requested terminal set.
   - `strict_road2edge`: one terminal per road node, with all road-2edge transformers/sources aggregated before OPT.
   - `contracted_corridor`: the older contracted ILP-corridor terminal skeleton, kept only for comparison.
2. Run balanced chain clustering and OPT structure construction for the requested `R`.
3. Export the OPT graph as a backbone solution with straight-line Euclidean edges.
4. Reuse the final-network builder to attach all remaining transformers as a street forest.
5. Run the deterministic switch-aware reliability decomposition with generalized chains.

## Implementation

- Script: `figures/optimal_hierarchy/opt_algorithm_sweep.py`
- Output table: `outputs\optimal_hierarchy\opt\p2u_opt_hierarchy_table.csv`
- Output JSON: `outputs\optimal_hierarchy\opt\p2u_opt_hierarchy_table.json`

Parameters:

```json
{
  "p_mean": 0.0005,
  "seed": 10,
  "kmeans_max_iter": 10,
  "strc_n_init_iters": 6,
  "strc_exact_vertices": true,
  "strc_trip_nearest_vertices": 24,
  "chain_n_init_iters": 6,
  "local_fix_max_changes": 0,
  "local_fix_max_risk_gain": 0.3,
  "generalized_method": "projection",
  "tree_mode": "street_forest",
  "target_source_contracted_r": true
}
```

Stage log:

```json
{
  "parameters": {
    "p_mean": 0.0005,
    "seed": 10,
    "kmeans_max_iter": 10,
    "strc_n_init_iters": 6,
    "strc_exact_vertices": true,
    "strc_trip_nearest_vertices": 24,
    "chain_n_init_iters": 6,
    "local_fix_max_changes": 0,
    "local_fix_max_risk_gain": 0.3,
    "generalized_method": "projection",
    "tree_mode": "street_forest",
    "target_source_contracted_r": true
  },
  "corridor_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_terminal_corridors_road2_k10_3857.gpkg",
  "terminal_summary": {
    "terminal_mode": "strict_road2edge",
    "road_graph_nodes": 71540,
    "road_graph_edges": 81206,
    "source_side_road_2edge_nodes": 58396,
    "transformers_total": 12373,
    "transformers_in_source_side_road_2edge": 10122,
    "sources_total": 15,
    "sources_in_source_side_road_2edge": 15,
    "aggregated_terminal_nodes": 9529,
    "aggregated_transformer_terminals": 9518,
    "aggregated_source_terminals": 15,
    "duplicate_road_nodes": 0,
    "capacity_kva": 879720.0
  },
  "backbone_terminal_count": 9529,
  "cases": {
    "50": {
      "backbone": {
        "algorithm": "OPT_equal_chain_backbone",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\opt\\R050\\p2u_opt_backbone_R50_3857.gpkg",
        "redundancy_constraint": 50,
        "physical_redundancy_constraint": 36,
        "target_source_contracted_r": true,
        "max_redundancy_constraint": null,
        "status_name": "OPT_CONSTRUCTED",
        "stop_reason": "opt_algorithm_completed",
        "runtime_s": 111.17252829996869,
        "seed": 10,
        "kmeans_max_iter": 10,
        "strc_n_init_iters": 6,
        "strc_exact_vertices": true,
        "strc_trip_nearest_vertices": 24,
        "chain_n_init_iters": 6,
        "local_fix": {
          "enabled": false
        },
        "input_nodes": 9529,
        "solution_nodes": 9529,
        "solution_edges": 9564,
        "solution_cycle_rank": 36,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "solution_transformer_nodes": 9514,
        "solution_source_nodes": 15,
        "input_duplicate_road_nodes": 0,
        "input_transformer_count": 10122,
        "input_source_count": 15,
        "objective_length_m": 522418.6922786822,
        "mean_edge_length_m": 54.62345172299061,
        "max_edge_length_m": 4213.923032388563,
        "tie_edges": 36,
        "tree_mode": "street_forest"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\opt\\R050\\p2u_opt_final_network_R50_risk.json",
      "case_runtime_s": 16.985095100011677
    },
    "100": {
      "backbone": {
        "algorithm": "OPT_equal_chain_backbone",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\opt\\R100\\p2u_opt_backbone_R100_3857.gpkg",
        "redundancy_constraint": 100,
        "physical_redundancy_constraint": 86,
        "target_source_contracted_r": true,
        "max_redundancy_constraint": null,
        "status_name": "OPT_CONSTRUCTED",
        "stop_reason": "opt_algorithm_completed",
        "runtime_s": 243.92847430007532,
        "seed": 10,
        "kmeans_max_iter": 10,
        "strc_n_init_iters": 6,
        "strc_exact_vertices": true,
        "strc_trip_nearest_vertices": 24,
        "chain_n_init_iters": 6,
        "local_fix": {
          "enabled": false
        },
        "input_nodes": 9529,
        "solution_nodes": 9529,
        "solution_edges": 9614,
        "solution_cycle_rank": 86,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "solution_transformer_nodes": 9514,
        "solution_source_nodes": 15,
        "input_duplicate_road_nodes": 0,
        "input_transformer_count": 10122,
        "input_source_count": 15,
        "objective_length_m": 575455.8826868907,
        "mean_edge_length_m": 59.85603106791041,
        "max_edge_length_m": 1291.2925557361693,
        "tie_edges": 86,
        "tree_mode": "street_forest"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\opt\\R100\\p2u_opt_final_network_R100_risk.json",
      "case_runtime_s": 18.947671300033107
    }
  },
  "reference_failure_rate_per_length": 1.0798531853768255e-05,
  "total_runtime_s": 99.90320319996681
}
```
