# P2U 3-Edge Projection Decomposition

## Results

**Goal**: Test whether the algebraic projection method from the sibling `cascading` project can replace the slow NetworkX 3-edge-component step for the full SMART-DS SFO P2U MV network.

**Main result**: The projection backend now runs through the full P2U network and produces generalized-chain data, but it is still expensive for repeated use: the measured decomposition runtime was `584.605` s after `GraphRel` initialization.

| metric | value |
| --- | ---: |
| Original MV nodes | 12572 |
| Original MV edges | 12728 |
| Sources | 15 |
| GraphRel nodes after source contraction | 12558 |
| GraphRel edges after source contraction | 12728 |
| Bridges | 8921 |
| 2-edge components | 8922 |
| Structure graph nodes | 1426 |
| Structure graph edges | 1590 |
| Regular chains | 1590 |
| 3-edge macro graph nodes | 1197 |
| 3-edge macro graph edges | 1421 |
| Generalized chains | 322 |
| Aggregated parallel macro edges | 61 |

The largest generalized chain lengths, in meters, were:

```text
4132.451, 2397.738, 2113.365, 2046.614, 1866.820,
1760.146, 1741.635, 1720.342, 1691.341, 1660.788
```

**Insights**: The projection method is mathematically usable on the P2U topology, and the real network does create generalized chains that are longer than the regular skeleton chains. However, the current implementation forms dense projection matrices on large analysis components, so it should be treated as an offline analysis backend unless we add more pruning or component-level caching.

**What this does not show**: This is not a reliability-risk decomposition yet. It only validates that the generalized-chain structure can be extracted from the full P2U MV topology. It also does not prove that projection is asymptotically faster for this distribution-network case.

**Reproduce**:

```powershell
conda run -n reliability python figures\p2u_three_edge_decomposition_benchmark.py --method projection
```

## Algorithm

The reliability graph first contracts all source nodes into the single `GraphRel.source` node. Bridges are removed to identify 2-edge-connected blocks. Inside the bridge-free part, maximal degree-2 chains are extracted as regular skeleton edges.

For generalized chains, the implementation builds an analysis graph from those regular chains. Each regular chain with endpoints `u,v` is represented as `u - a_i - v`, where `a_i = ("reduced_edge", i)` is an auxiliary node. This preserves parallel reduced paths while keeping a simple `nx.Graph`.

The projection backend calls the existing cascading detector:

```python
dc_graph.structure._raw_components_projection(analysis_graph, k=3, min_skeleton_nodes=2)
```

Conceptually, this computes 3-edge-connected skeleton components by detecting near-zero two-edge cuts using the edge projection matrix

$$
P_{ef} = \sqrt{b_e}\, h_e^T L^+ h_f\, \sqrt{b_f}.
$$

After 3-edge cores are found, each core is contracted to a macro node. Repeated macro edges are aggregated into one simple edge with summed `length`, summed `edge_weight`, provenance metadata, and `parallel_macro_edge_count`. This keeps the user-facing decomposition in simple-graph form while preserving the fact that several reduced chains may connect the same macro blocks.

## Implementation

The public API is:

```python
GraphRel.decompose(
    include_generalized_chains=True,
    generalized_component_method="projection",
)
```

Relevant files:

- `indexes/graph_rel.py`: exposes the `generalized_component_method` option.
- `indexes/rel_decomposition.py`: builds the chain-based analysis graph, loads the cascading projection detector, contracts 3-edge components, and aggregates parallel macro edges.
- `figures/p2u_three_edge_decomposition_benchmark.py`: reproducible P2U benchmark entry point.
- `tests/test_rel_decomposition.py`: validates that the projection backend matches NetworkX on a generated graph with three dense blobs connected by two-cut structure.

The implementation still rejects user-supplied `nx.MultiGraph` inputs. The only parallel handling added here is internal aggregation after contracting 3-edge macro components, because the P2U simple physical graph naturally induces repeated macro links.

Known risks:

- The projection detector is imported from the sibling `cascading` project through its internal function. This is intentionally minimal, but it means API changes in `cascading` can break the optional backend.
- The full P2U projection run is slow because dense projection blocks are large.
- No project research log was found at `codex/documents/research_logs/research_log.md`, so there was no research-log entry to update.
