# Tie-Switching Simulation Plan

Date: 2026-08-05
Branch: `codex/tie-switching-simulations`

## Goal

Add switch-aware reliability simulations for distribution networks, especially SMART-DS power networks, without changing the core exact reliability/tree-decomposition algorithms. The manuscript currently reports MCMC reliability simulations, so the implementation should focus on the MCMC path.

## Current Code Structure

Relevant files:

- `indexes/simulation.py`: event-driven MCMC reliability simulator.
- `indexes/graph_rel.py`: `GraphRel.calc_rel_simulation()` wraps `simulate_rel()`.
- `figures/real_data.py`: builds real communication/water/power networks and computes original-vs-optimal reliability.
- `utilities/figures_utilities.py`: helper `saidi_with_lengths()` used by several simulation/figure scripts.
- `optimal_network/other_methods/improve_tree.py`: greedy tree-improvement method that adds chords; these chords can be interpreted as tie switches.
- `optimal_network/`: exact/constructive optimal-network machinery. This does not need to change for the switch-aware MCMC task.

Important current behavior:

- `figures/real_data.py::get_sfo_bounded_network()` builds a SMART-DS graph, then immediately reduces it to a minimum spanning tree:

```python
G = nx.edge_subgraph(G, nx.minimum_spanning_tree(G).edges).copy()
```

- `figures/real_data.py::get_sfo_networks(improve=True)` then calls `improve_tree()` to add redundant edges.
- `figures/real_data.py::get_G_and_OG_saidi()` computes length-based edge failure probabilities and calls `GraphRel.calc_rel_simulation()`.
- `indexes/simulation.py::ConnectedComponents.reliability()` currently evaluates reliability by connected component in the live graph after failed edges are removed.

This means the existing MCMC already answers: "which loads are connected in the energized graph after edge failures?" It does not yet distinguish normal radial operation from a larger available restoration graph with normally open tie switches.

## Design Decision

Do not create a completely separate simulation framework.

Instead, add a small switch-aware layer to the existing MCMC simulator:

- Keep `simulate_rel()` as the core event engine.
- Add an optional restoration-policy object/function used by `ConnectedComponents.reliability()`.
- Keep normal behavior as the default.
- Add separate real-data functions/figure scripts for the switch-aware Nature Energy validation.

This keeps the old figures reproducible while allowing switch-aware reliability where needed.

## Theory-To-Code Mapping

Theory graph definitions:

- `T`: normally closed radial graph.
- `E_tie`: normally open redundant/tie edges.
- `G_avail = T union E_tie`: physical graph available after switching.
- `X`: failed edges.
- `G_live = G_avail \ X`: graph after failures.
- A load is served if it can be connected to a source in `G_live` after allowed switching.

Under ideal fast switching and no power-flow constraints, restoration is a graph connectivity/radialization problem.

## Algorithms Needed

### 1. Fast Connectivity Restoration

Purpose:

Compute the best possible post-fault served load if switching is ideal and electrical constraints are ignored.

Algorithm:

1. Remove failed edges from the available graph.
2. Find all nodes connected to at least one source.
3. Served load is the total weight in source-connected components.
4. The radial operating graph can be chosen as any spanning forest of the source-connected part.

This is enough for the main theory because if a node is connected in `G_avail \ X`, a radial path to a source exists inside a spanning tree of that component.

Expected function:

```python
restored_connected_weight(G_available, failed_edges, sources, weight_attr="weight")
```

Complexity:

- `O(N + M)` per event-state reliability recomputation.
- Compatible with the current MCMC design.

### 2. Optimal Real-Time Radial Restoration

Purpose:

Choose a radial switching configuration after faults that maximizes served weighted load.

For the pure connectivity model, this is simple:

1. Build `G_live = G_available \ failed_edges`.
2. Add a super-source connected to all real sources.
3. Keep all components connected to the super-source.
4. Return a spanning tree/forest of the served subgraph.

If all served loads have no electrical constraints, every connected load can be served, so the objective is solved by connectivity alone.

Expected function:

```python
optimal_radial_restoration(
    G_available,
    failed_edges,
    sources,
    weight_attr="weight",
    edge_cost_attr=None,
) -> tuple[nx.Graph, set]
```

Return:

- restored radial graph,
- served nodes or disconnected nodes.

### 3. Optional Constrained Restoration

Purpose:

Prepare a later SI extension if reviewers demand realism beyond connectivity.

Possible constraints:

- one source per island,
- no cycles in energized state,
- no switching of non-switchable lines,
- voltage/ampacity constraints checked by OpenDSS,
- maximum number of switching actions,
- priority load weights.

This should not be in the first implementation unless needed. The first Nature Energy figure can use ideal fast switching and state the assumption.

### 4. Switch-Aware Effective Risk

Purpose:

Estimate the switch-aware `\tilde{\lambda}_{sw}` used in the SI theory.

Algorithm:

1. For each region/chain or for the whole graph, enumerate sampled two-edge failures or all two-edge failures if small enough.
2. For each pair, compute post-restoration disconnected load.
3. Sum `lambda_i lambda_j w_ij`.

Expected function:

```python
estimate_switch_effective_risk(
    G_available,
    sources,
    edge_failure_weights,
    weight_attr="weight",
    max_pairs=None,
    rng=None,
)
```

This can support a small SI validation figure.

## Data Model Changes

Graph-level attributes:

- `graph["sources"]`: existing source list.
- `graph["normal_edges"]`: edges normally closed in radial operation.
- `graph["tie_edges"]`: normally open edges that can be closed after faults.
- `graph["available_graph_kind"]`: optional label such as `"smart_ds_switches"` or `"tree_plus_added_ties"`.

Edge attributes:

- `is_switch`: whether the line/device is switchable.
- `is_open`: whether it is normally open.
- `normally_closed`: whether energized in normal operation.
- `is_tie`: whether it is a normally open restoration edge.
- `length_m` or `length`: used for failure probability.
- `prob`: current failure probability used by MCMC.

For optimized networks:

- Original tree edges: `normally_closed=True`, `is_tie=False`.
- Added redundant edges: `normally_closed=False`, `is_tie=True`, `is_switch=True`.

## SMART-DS Loader Changes

Current loader loses switch information because it uses only `Line_N.shp` and then MST-reduces the graph.

Needed changes:

1. Add a loader that builds both:
   - normal radial graph,
   - available graph with tie/switch edges.

2. Use `SwitchingDevices_N.shp` where possible:
   - fields seen from file header include `NodeA`, `NodeB`, `Code`, `NomV_kV`, `Subest`, `Feeder`;
   - records include breakers and pad switches.

3. Use `data.json`/GeoJSON properties if needed:
   - observed properties include `is_switch`, `is_open`, `is_fuse`, `length (km)`, `ampacity`, `phases`, `type`, `name`.

4. Preserve the current tree-only graph for baseline/original comparisons.

Possible new function names:

```python
get_sfo_switching_network(area: str, box: str) -> tuple[nx.Graph, nx.Graph]
build_mv_switching_network_from_boundary(data_dir, boundary)
```

## MCMC Integration Plan

Add an optional parameter to `GraphRel.calc_rel_simulation()` and `simulate_rel()`:

```python
restoration_policy: RestorationPolicy | None = None
```

Where `RestorationPolicy` can be a callable:

```python
Callable[[nx.Graph, set[tuple], object, str, RelType], float]
```

For cleaner code, better define a small class in a new file:

```python
indexes/restoration.py
```

Suggested contents:

- `NoRestorationPolicy`
- `ConnectivityRestorationPolicy`
- `RadialRestorationResult`
- `optimal_radial_restoration()`
- `restored_disconnected_fraction()`

Then `ConnectedComponents.reliability()` can delegate to the policy when provided.

## Figure/Data Plan

Create one new real-data simulation output, rather than modifying every old figure.

Suggested new output:

- `data/real_networks_switching.nxjson`
- `outputs/real_data_switching/real_data_switching.svg`

Suggested comparisons:

1. Original radial feeder, no restoration.
2. Original feeder plus existing SMART-DS switches/ties, switch-aware restoration.
3. Optimized tree plus added tie switches, switch-aware restoration.

Metrics:

- `saidi_no_switching`
- `saidi_existing_switching`
- `saidi_optimized_switching`
- `rel_gain_existing_switching`
- `rel_gain_optimized_switching`
- total network length/cost
- redundancy count `R`
- switch/tie count

## What Not To Change

Do not edit the exact reliability/tree-decomposition machinery for this phase:

- `indexes/tree_decomposition.py`
- `indexes/partition_collection.py`
- exact polynomial reliability code in `indexes/graph_rel.py`
- core optimal network construction in `optimal_network/construct_strc.py`, `cluster_chains.py`, `connect_chains.py`

Those support the theory and existing figures. The switch-aware result can be added through MCMC and real-data validation.

## Recommended Implementation Order

1. Add `indexes/restoration.py` with pure graph restoration functions and small tests/smoke examples.
2. Extend `indexes/simulation.py` so `ConnectedComponents` can evaluate reliability using a restoration policy.
3. Extend `GraphRel.calc_rel_simulation()` to pass the policy through.
4. Add a helper in `utilities/figures_utilities.py`, probably `saidi_with_lengths_and_switching()`.
5. Add SMART-DS switch-aware graph construction in `figures/real_data.py`.
6. Add a new figure script or new mode in `figures/real_data_plot.py`.
7. Run a small synthetic sanity check:
   - tree line: single edge failure disconnects downstream load,
   - tree plus tie: one internal failure can be restored,
   - two failures can still disconnect internal blocks.
8. Run the SFO Pacific/Davidson switch-aware simulation.

## Expected Manuscript Claim After Coding

After this implementation, the paper can honestly say:

> For radially operated distribution networks, we evaluate sustained interruption after ideal post-fault switching. Redundant lines are treated as normally open tie switches; after each failure event, the model computes the maximum load that can be reconnected to a source by radial reconfiguration. The resulting switch-aware reliability preserves the predicted cost-reliability scaling while reducing the effective internal chain risk.

