# P2U Hierarchical Redundancy Sweep (exact)

## Results

**Goal**: prepare an article-level P2U experiment where a road-constrained 2-edge-connected backbone is built under a redundancy target and all remaining transformers are attached as a street forest.

**Main result**: this manifest is the first clean sweep scaffold. It uses the R50-style parameters that produced the high `O(p^2)` structural-risk example and records the achieved network for each requested redundancy budget.

All reliability rows use deterministic decomposition with length-scaled `p_mean = 5e-4`.

| network | r_request | length_km | r_theory | z_w | z_r | z_f | z_f_p | z_f_p2 | risk_total | risk_o_p | risk_o_p2 | risk_tree | risk_section | risk_internal | risk_structural | bridges | generalized_chains | gen_lambda_mean_km | gen_lambda_sigma_over_mean | gen_lambda_max_km | ilp_runtime_s | decomposition_runtime_s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Original P2U MV |  | 5.89e2 | 171 | 0 | 0 | 0 | 0 | 0 | 2.27e-3 | 2.26e-3 | 6.15e-6 | 2.26e-3 | 0 | 3.61e-7 | 5.79e-6 | 8921 | 324 | 3.24e-1 | 1.45e0 | 4.04e0 |  | 7.72e0 |
| Hierarchical road R50 | 50 | 5.88e2 | 50 | 2.63e-3 | 7.08e-1 | 7.03e-1 | 8.19e-1 | -4.21e1 | 6.74e-4 | 4.09e-4 | 2.65e-4 | 3.41e-4 | 6.76e-5 | 6.39e-6 | 2.59e-4 | 5206 | 123 | 2.93e0 | 1.29e0 | 1.84e1 | 5.73e2 | 5.60e0 |
| Hierarchical road R70 | 70 | 5.88e2 | 70 | 2.09e-3 | 5.91e-1 | 7.32e-1 | 8.19e-1 | -3.12e1 | 6.08e-4 | 4.09e-4 | 1.98e-4 | 3.42e-4 | 6.77e-5 | 6.39e-6 | 1.92e-4 | 5206 | 148 | 2.48e0 | 1.38e0 | 1.77e1 | 1.08e0 | 6.17e0 |
| Hierarchical road R90 | 90 | 5.89e2 | 90 | 1.32e-3 | 4.74e-1 | 7.34e-1 | 8.19e-1 | -3.03e1 | 6.02e-4 | 4.10e-4 | 1.92e-4 | 3.42e-4 | 6.78e-5 | 6.40e-6 | 1.86e-4 | 5206 | 172 | 2.25e0 | 1.56e0 | 1.90e1 | 5.52e-1 | 6.71e0 |
| Hierarchical road R110 | 110 | 5.89e2 | 110 | 4.10e-4 | 3.57e-1 | 7.41e-1 | 8.18e-1 | -2.78e1 | 5.88e-4 | 4.10e-4 | 1.77e-4 | 3.42e-4 | 6.79e-5 | 6.34e-6 | 1.71e-4 | 5206 | 202 | 2.00e0 | 1.74e0 | 2.06e1 | 6.66e-1 | 8.99e0 |
| Hierarchical road R130 | 130 | 5.90e2 | 130 | -5.62e-4 | 2.40e-1 | 7.50e-1 | 8.18e-1 | -2.45e1 | 5.67e-4 | 4.11e-4 | 1.57e-4 | 3.43e-4 | 6.79e-5 | 6.33e-6 | 1.51e-4 | 5206 | 232 | 1.78e0 | 1.88e0 | 2.21e1 | 9.74e-1 | 8.90e0 |
| Hierarchical road R150 | 150 | 5.88e2 | 150 | 1.71e-3 | 1.23e-1 | 7.90e-1 | 8.18e-1 | -9.37e0 | 4.76e-4 | 4.12e-4 | 6.38e-5 | 3.45e-4 | 6.66e-5 | 5.98e-6 | 5.79e-5 | 5224 | 332 | 1.09e0 | 1.91e0 | 1.62e1 | 6.56e0 | 8.16e0 |
| Hierarchical road R170 | 170 | 5.89e2 | 170 | 9.35e-4 | 5.85e-3 | 7.90e-1 | 8.18e-1 | -9.18e0 | 4.75e-4 | 4.12e-4 | 6.26e-5 | 3.46e-4 | 6.67e-5 | 5.98e-6 | 5.67e-5 | 5224 | 358 | 1.04e0 | 2.02e0 | 1.68e1 | 1.17e0 | 9.53e0 |
| Hierarchical road R171 | 171 | 5.89e2 | 171 | 1.17e-3 | 0 | 7.92e-1 | 8.18e-1 | -8.50e0 | 4.71e-4 | 4.12e-4 | 5.85e-5 | 3.45e-4 | 6.69e-5 | 5.97e-6 | 5.25e-5 | 5221 | 373 | 9.93e-1 | 2.05e0 | 1.72e1 | 5.64e0 | 1.01e1 |

**Figures**:

- `Hierarchical road R50`: ![](R050/p2u_final_network_R50_map.png)
- `Hierarchical road R70`: ![](R070/p2u_final_network_R70_map.png)
- `Hierarchical road R90`: ![](R090/p2u_final_network_R90_map.png)
- `Hierarchical road R110`: ![](R110/p2u_final_network_R110_map.png)
- `Hierarchical road R130`: ![](R130/p2u_final_network_R130_map.png)
- `Hierarchical road R150`: ![](R150/p2u_final_network_R150_map.png)
- `Hierarchical road R170`: ![](R170/p2u_final_network_R170_map.png)
- `Hierarchical road R171`: ![](R171/p2u_final_network_R171_map.png)

**Insights**:

- `O(p) = Tree + Section` is the first-order bridge/section risk.
- `O(p^2) = Internal + Structural` is the second-order chain and generalized-chain risk.
- `gen_lambda_*` reports demand-normalized generalized-chain effective lengths, using `tilde_lambda_q = lambda_q sqrt(Q w_q / W)`.
- The R50 case is kept as the first benchmark because it shows that a low-redundancy 2-connected backbone can leave a large `O(p^2)` component.
- This run enforces exact source-contracted redundancy, so `r_request` should match `r_theory` for successful rows.


**What this does not show**:

- This stage produces a simulation table and preview maps, not the final article figure layout.
- ILP results are only as strong as the recorded Gurobi status and MIP gap for each row.
- The experiment remains connectivity reliability, not voltage/power-flow validation.

**Reproduce**:

```powershell
& C:\Users\rotem\anaconda3\envs\reliability\python.exe figures\optimal_hierarchy\redundancy_sweep.py --redundancy-mode exact --r-values 50 70 90 110 130 150 170 171 --reuse-existing
```

## Algorithm

For each requested redundancy value, the pipeline uses the road-corridor terminal graph and solves a minimum-length 2-edge-connected backbone ILP:

$$
\min \sum_e w_e x_e
$$

subject to degree, source-incidence, lazy 2-edge cut constraints, and the selected redundancy constraint:

$$
|E_\mathrm{selected}| - |V| + 1 = R_\mathrm{target}.
$$

The final graph attaches non-backbone transformer terminals using shortest paths on the street network and contracts street/transformer chains into switch-section edges.

The reliability split is:

$$
F = F_{O(p)} + F_{O(p^2)}, \quad F_{O(p)} = F_\mathrm{tree} + F_\mathrm{section}, \quad F_{O(p^2)} = F_\mathrm{internal} + F_\mathrm{structural}.
$$

Relative indexes are computed against the original P2U MV network:

$$
Z_W = 1 - W/W_0, \quad Z_R = 1 - R/R_0, \quad Z_F = 1 - F/F_0.
$$

Current parameters:

- `p_mean = 0.0005`
- `time_limit = 600.0` s
- `MIPGap = 0.05`
- `Threads = 0`
- `cut_mode = callback`
- `tree_mode = street_forest`
- `generalized_method = projection`
- `redundancy_mode = exact`

## Implementation

Entry point:

- `figures/optimal_hierarchy/redundancy_sweep.py`

Shared code reused from the exploratory P2U implementation:

- `figures/optimal_sfo/prepare_p2u_corridor_network.py` builds the road-corridor candidate graph.
- `figures/optimal_sfo/run_p2u_ilp_2edge.py` solves the backbone ILP.
- `figures/optimal_sfo/p2u_final_network.py` builds the final backbone-plus-forest graph.
- `figures/optimal_sfo/analyze_p2u_final_network.py` computes topology, load, and length summaries.
- `figures/optimal_sfo/compare_p2u_old_new_reliability.py` provides the deterministic risk-decomposition path.

Outputs:

- `outputs\optimal_hierarchy\exact\p2u_hierarchical_redundancy_sweep_exact_table.csv`
- `outputs\optimal_hierarchy\exact\p2u_hierarchical_redundancy_sweep_exact_table.json`
- `outputs\optimal_hierarchy\exact\p2u_hierarchical_redundancy_sweep_exact.md`
- `outputs\optimal_hierarchy\exact\p2u_hierarchical_redundancy_sweep_exact_diagnostic.png`
- `outputs\optimal_hierarchy\exact/R*/p2u_final_network_*_map.png`

Stage log:

```json
{
  "parameters": {
    "p_mean": 0.0005,
    "ilp_time_limit_s": 600.0,
    "ilp_mip_gap": 0.05,
    "ilp_threads": 0,
    "ilp_max_cut_rounds": 100,
    "ilp_cut_mode": "callback",
    "tree_mode": "street_forest",
    "generalized_method": "projection",
    "redundancy_mode": "exact"
  },
  "cases": {
    "50": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 572.9616575241089,
        "cut_rounds": 144,
        "added_cuts": 13491,
        "objective_length_m": 400041.9126132575,
        "best_bound": 396472.2314545691,
        "mip_gap": 0.008923267903029431,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5345,
        "solution_cycle_rank": 50,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": null,
        "max_redundancy_constraint": 50,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R050\\p2u_backbone_R50_3857.gpkg",
        "imported_from": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_ilp_2edge_solution_Rmax50_summary.json",
        "redundancy_mode": "exact",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R050\\p2u_backbone_R50_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6616,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 187749.2901489289,
        "street_branch_nodes": 267,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R050\\p2u_final_network_R50_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R050\\p2u_final_network_R50_streetforest_metadata.json",
        "imported_from": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_final_network_Rmax50_streetforest_metadata.json",
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R050\\p2u_final_network_R50_risk.json",
      "plot": null,
      "case_runtime_s": 10.020475900033489
    },
    "70": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 1.083369255065918,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 400361.1005312399,
        "best_bound": 380401.2103711625,
        "mip_gap": 0.04985471898641639,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5365,
        "solution_cycle_rank": 70,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 70,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5345,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R050\\p2u_backbone_R50_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R070\\p2u_backbone_R70_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R070\\p2u_backbone_R70_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6616,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 187749.2901489289,
        "street_branch_nodes": 267,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R070\\p2u_final_network_R70_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R070\\p2u_final_network_R70_streetforest_metadata.json",
        "stage_runtime_s": 71.034295200021,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R070\\p2u_final_network_R70_risk.json",
      "plot": null,
      "case_runtime_s": 10.787777799996547
    },
    "90": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 0.552309513092041,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 400812.2718919704,
        "best_bound": 380772.2016643274,
        "mip_gap": 0.04999864433553158,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5385,
        "solution_cycle_rank": 90,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 90,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5365,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R070\\p2u_backbone_R70_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R090\\p2u_backbone_R90_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R090\\p2u_backbone_R90_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6616,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 187749.2901489289,
        "street_branch_nodes": 267,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R090\\p2u_final_network_R90_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R090\\p2u_final_network_R90_streetforest_metadata.json",
        "stage_runtime_s": 34.99637750000693,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R090\\p2u_final_network_R90_risk.json",
      "plot": null,
      "case_runtime_s": 10.895551900030114
    },
    "110": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 0.6660184860229492,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 401348.56797050894,
        "best_bound": 381319.5278737027,
        "mip_gap": 0.04990435171623181,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5405,
        "solution_cycle_rank": 110,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 110,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5385,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R090\\p2u_backbone_R90_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R110\\p2u_backbone_R110_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R110\\p2u_backbone_R110_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6616,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 187749.2901489289,
        "street_branch_nodes": 267,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R110\\p2u_final_network_R110_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R110\\p2u_final_network_R110_streetforest_metadata.json",
        "stage_runtime_s": 40.23630150000099,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R110\\p2u_final_network_R110_risk.json",
      "plot": null,
      "case_runtime_s": 13.3920378999901
    },
    "130": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 0.9743261337280273,
        "cut_rounds": 1,
        "added_cuts": 19,
        "objective_length_m": 401921.1090319074,
        "best_bound": 381843.8049967757,
        "mip_gap": 0.04995334552965174,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5425,
        "solution_cycle_rank": 130,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 130,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5405,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R110\\p2u_backbone_R110_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R130\\p2u_backbone_R130_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R130\\p2u_backbone_R130_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6616,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 187749.2901489289,
        "street_branch_nodes": 267,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R130\\p2u_final_network_R130_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R130\\p2u_final_network_R130_streetforest_metadata.json",
        "stage_runtime_s": 35.280984800017904,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R130\\p2u_final_network_R130_risk.json",
      "plot": null,
      "case_runtime_s": 13.371031999995466
    },
    "150": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 6.561850309371948,
        "cut_rounds": 12,
        "added_cuts": 4442,
        "objective_length_m": 399974.4463413936,
        "best_bound": 382189.28330135165,
        "mip_gap": 0.04446574825648137,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5445,
        "solution_cycle_rank": 150,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 150,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5395,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R100\\p2u_backbone_R100_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_backbone_R150_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_backbone_R150_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6602,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 188359.00744761626,
        "street_branch_nodes": 271,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_final_network_R150_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_final_network_R150_streetforest_metadata.json",
        "stage_runtime_s": 20.384308099979535,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_final_network_R150_risk.json",
      "plot": null,
      "case_runtime_s": 11.673828699975275
    },
    "170": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 1.165802001953125,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 400429.12416411674,
        "best_bound": 380410.47662493534,
        "mip_gap": 0.0499929858522896,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5465,
        "solution_cycle_rank": 170,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 170,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5445,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_backbone_R150_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R170\\p2u_backbone_R170_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R170\\p2u_backbone_R170_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6602,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 188359.00744761626,
        "street_branch_nodes": 271,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R170\\p2u_final_network_R170_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R170\\p2u_final_network_R170_streetforest_metadata.json",
        "stage_runtime_s": 70.1345073999837,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R170\\p2u_final_network_R170_risk.json",
      "plot": null,
      "case_runtime_s": 14.190810200001579
    },
    "171": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 5.642043590545654,
        "cut_rounds": 15,
        "added_cuts": 3735,
        "objective_length_m": 400349.1246976708,
        "best_bound": 382642.39927402185,
        "mip_gap": 0.04422821065743664,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5466,
        "solution_cycle_rank": 171,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": 171,
        "max_redundancy_constraint": null,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5445,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R150\\p2u_backbone_R150_3857.gpkg",
        "redundancy_mode": "exact",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R171\\p2u_backbone_R171_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R171\\p2u_backbone_R171_3857.gpkg",
        "final_graph_definition": "switch-section graph: graph nodes are transformer/source terminals plus optional zero-load street branch nodes; degree-2 transformer chains are stored as edge load; tree mode is street_forest",
        "tree_mode": "street_forest",
        "retained_backbone_road_nodes": 5310,
        "all_backbone_terminal_sequence_road_nodes": 6605,
        "total_transformers_in_data": 12373,
        "unattached_tree_terminal_count": 0,
        "tree_attachment_physical_road_union_length_m": 188303.66548790035,
        "street_branch_nodes": 271,
        "lv_assignment": {
          "lv_components": 12373,
          "lv_nodes": 69004,
          "lv_edges": 56631,
          "lv_components_without_transformer": 0,
          "lv_components_with_multiple_transformers": 0,
          "consumer_points_total": 46379,
          "consumer_points_assigned_to_mv": 46379,
          "consumer_points_unassigned": 0,
          "consumer_points_assigned_by_nearest_transformer_fallback": 196,
          "nearest_transformer_fallback_mean_distance_m": 44.595398456044656,
          "nearest_transformer_fallback_max_distance_m": 203.9177052137165,
          "demand_kw_total_raw": 661466.7399999977,
          "demand_kw_assigned_to_mv": 661466.74,
          "num_customers_total_raw": 191140.0,
          "num_customers_assigned_to_mv": 191140.0,
          "yearly_kwh_total_raw": 5876669360.610104,
          "yearly_kwh_assigned_to_mv": 5876669360.61,
          "transformer_nodes_with_assigned_consumers": 12373
        },
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R171\\p2u_final_network_R171_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R171\\p2u_final_network_R171_streetforest_metadata.json",
        "stage_runtime_s": 22.07461080001667,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\exact\\R171\\p2u_final_network_R171_risk.json",
      "plot": null,
      "case_runtime_s": 15.382748599979095
    }
  },
  "corridors": {
    "status": "reused",
    "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_terminal_corridors_road2_k10_3857.gpkg"
  },
  "total_runtime_s": 114.95644789998187
}
```

No project research log was found at `codex/documents/research_logs/research_log.md`, so no research-log entry was updated.
