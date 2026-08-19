# P2U Hierarchical Redundancy Sweep

## Results

**Goal**: prepare an article-level P2U experiment where a road-constrained 2-edge-connected backbone is built under a redundancy budget and all remaining transformers are attached as a street forest.

**Main result**: this manifest is the first clean sweep scaffold. It uses the R50-style parameters that produced the high `O(p^2)` structural-risk example and records the achieved network for each requested redundancy budget.

All reliability rows use deterministic decomposition with length-scaled `p_mean = 5e-4`.

| network | r_max | length_km | r_theory | z_w | z_r | z_f | risk_total | risk_o_p | risk_o_p2 | risk_tree | risk_section | risk_internal | risk_structural | bridges | ilp_runtime_s | decomposition_runtime_s |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Original P2U MV |  | 589.339 | 171 | 0 | 0 | 0 | 0.0022664 | 0.00226025 | 6.15404e-06 | 0.00226025 | 0 | 3.6127e-07 | 5.79277e-06 | 8921 |  | 7.07645 |
| Hierarchical road R50 | 50 | 587.791 | 50 | 0.00262696 | 0.707602 | 0.702528 | 0.000674192 | 0.000408827 | 0.000265365 | 0.000341193 | 6.76344e-05 | 6.38796e-06 | 0.000258977 | 5206 | 572.962 | 5.49015 |
| Hierarchical road R100 | 100 | 587.791 | 50 | 0.00262696 | 0.707602 | 0.702525 | 0.000674198 | 0.000408827 | 0.000265371 | 0.000341193 | 6.76344e-05 | 6.3875e-06 | 0.000258984 | 5206 | 0.372157 | 5.40761 |
| Hierarchical road R150 | 150 | 587.791 | 50 | 0.00262696 | 0.707602 | 0.702525 | 0.000674198 | 0.000408827 | 0.000265371 | 0.000341193 | 6.76344e-05 | 6.3875e-06 | 0.000258984 | 5206 | 0.229868 | 14.4853 |
| Hierarchical road R171 | 171 | 587.791 | 50 | 0.00262696 | 0.707602 | 0.702525 | 0.000674198 | 0.000408827 | 0.000265371 | 0.000341193 | 6.76344e-05 | 6.3875e-06 | 0.000258984 | 5206 | 0.358159 | 13.7327 |

**Figures**:

- `Hierarchical road R50`: ![](R050/p2u_final_network_Rmax50_map.png)
- `Hierarchical road R100`: ![](R100/p2u_final_network_Rmax100_map.png)
- `Hierarchical road R150`: ![](R150/p2u_final_network_Rmax150_map.png)
- `Hierarchical road R171`: ![](R171/p2u_final_network_Rmax171_map.png)

**Insights**:

- `O(p) = Tree + Section` is the first-order bridge/section risk.
- `O(p^2) = Internal + Structural` is the second-order chain and generalized-chain risk.
- The R50 case is kept as the first benchmark because it shows that a low-redundancy 2-connected backbone can leave a large `O(p^2)` component.
- All requested upper bounds selected the same achieved `R=50` backbone. This is expected for a pure minimum-length objective with `R <= R_max`: extra redundancy is allowed but not rewarded.

**What this does not show**:

- This stage produces a simulation table and preview maps, not the final article figure layout.
- ILP results are only as strong as the recorded Gurobi status and MIP gap for each row.
- The experiment remains connectivity reliability, not voltage/power-flow validation.

**Reproduce**:

```powershell
& C:\Users\rotem\anaconda3\envs\reliability\python.exe figures\optimal_hierarchy\redundancy_sweep.py --r-values 50 100 150 171 --reuse-existing
```

## Algorithm

For each redundancy budget `R_max`, the pipeline uses the road-corridor terminal graph and solves a minimum-length 2-edge-connected backbone ILP:

$$
\min \sum_e w_e x_e
$$

subject to degree, source-incidence, lazy 2-edge cut constraints, and:

$$
|E_\mathrm{selected}| - |V| + 1 \le R_\max.
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

- `outputs/optimal_hierarchy/p2u_hierarchical_redundancy_sweep_table.csv`
- `outputs/optimal_hierarchy/p2u_hierarchical_redundancy_sweep_table.json`
- `outputs/optimal_hierarchy/p2u_hierarchical_redundancy_sweep.md`
- `outputs/optimal_hierarchy/R*/p2u_final_network_Rmax*_map.png`

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
    "generalized_method": "projection"
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
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_backbone_Rmax50_3857.gpkg",
        "imported_from": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_ilp_2edge_solution_Rmax50_summary.json",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_backbone_Rmax50_3857.gpkg",
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
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_final_network_Rmax50_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_final_network_Rmax50_streetforest_metadata.json",
        "imported_from": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_final_network_Rmax50_streetforest_metadata.json",
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_final_network_Rmax50_risk.json",
      "plot": {
        "stage_status": "reused",
        "output_png": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_final_network_Rmax50_map.png"
      },
      "case_runtime_s": 0.0006150999688543379
    },
    "100": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 0.37215733528137207,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 400041.9126132575,
        "best_bound": 380432.13069042185,
        "mip_gap": 0.049019318487744254,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5345,
        "solution_cycle_rank": 50,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": null,
        "max_redundancy_constraint": 100,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5345,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R050\\p2u_backbone_Rmax50_3857.gpkg",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_backbone_Rmax100_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_backbone_Rmax100_3857.gpkg",
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
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_final_network_Rmax100_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_final_network_Rmax100_streetforest_metadata.json",
        "stage_runtime_s": 35.05665789998602,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_final_network_Rmax100_risk.json",
      "plot": {
        "stage_status": "reused",
        "output_png": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_final_network_Rmax100_map.png"
      },
      "case_runtime_s": 0.0005490999901667237
    },
    "150": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 0.22986793518066406,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 400041.9126132575,
        "best_bound": 380432.13069042185,
        "mip_gap": 0.049019318487744254,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5345,
        "solution_cycle_rank": 50,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": null,
        "max_redundancy_constraint": 150,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5345,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R100\\p2u_backbone_Rmax100_3857.gpkg",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_backbone_Rmax150_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_backbone_Rmax150_3857.gpkg",
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
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_final_network_Rmax150_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_final_network_Rmax150_streetforest_metadata.json",
        "stage_runtime_s": 24.487962500017602,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_final_network_Rmax150_risk.json",
      "plot": {
        "stage_status": "reused",
        "output_png": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_final_network_Rmax150_map.png"
      },
      "case_runtime_s": 0.0005378000205382705
    },
    "171": {
      "backbone": {
        "status": 2,
        "status_name": "OPTIMAL",
        "stop_reason": "2edge_connected",
        "runtime_s": 0.35815930366516113,
        "cut_rounds": 0,
        "added_cuts": 0,
        "objective_length_m": 400041.9126132575,
        "best_bound": 380432.13069042185,
        "mip_gap": 0.049019318487744254,
        "input_nodes": 5296,
        "input_edges": 16985,
        "solution_nodes": 5296,
        "solution_edges": 5345,
        "solution_cycle_rank": 50,
        "solution_connected": true,
        "solution_bridge_count": 0,
        "solution_is_2edge_connected": true,
        "redundancy_constraint": null,
        "max_redundancy_constraint": 171,
        "time_limit_s": 600.0,
        "threads": 0,
        "requested_mip_gap": 0.05,
        "cut_mode": "callback",
        "source_connectivity_constraints": 15,
        "warm_start_edges": 5345,
        "selected_original_sources_with_incident_edge": 15,
        "original_sources_with_candidate_edges": 15,
        "warm_start_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R150\\p2u_backbone_Rmax150_3857.gpkg",
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R171\\p2u_backbone_Rmax171_3857.gpkg",
        "stage_status": "reused"
      },
      "final_network": {
        "input_backbone_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R171\\p2u_backbone_Rmax171_3857.gpkg",
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
        "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R171\\p2u_final_network_Rmax171_streetforest_3857.gpkg",
        "metadata_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R171\\p2u_final_network_Rmax171_streetforest_metadata.json",
        "stage_runtime_s": 52.760943299974315,
        "stage_status": "reused"
      },
      "risk_json": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R171\\p2u_final_network_Rmax171_risk.json",
      "plot": {
        "stage_status": "reused",
        "output_png": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_hierarchy\\R171\\p2u_final_network_Rmax171_map.png"
      },
      "case_runtime_s": 0.0005301000201143324
    }
  },
  "corridors": {
    "status": "reused",
    "output_gpkg": "C:\\Users\\rotem\\Desktop\\\u05de\u05e1\u05de\u05db\u05d9\u05dd\\\u05ea\u05d5\u05d0\u05e8\\\u05ea\u05d6\u05d4\\\u05e1\u05d9\u05de\u05d5\u05dc\u05e6\u05d9\u05d5\u05ea\\cost_effective\\outputs\\optimal_sfo\\p2u_terminal_corridors_road2_k10_3857.gpkg"
  },
  "total_runtime_s": 0.0030200000037439167
}
```

No project research log was found at `codex/documents/research_logs/research_log.md`, so no research-log entry was updated.
