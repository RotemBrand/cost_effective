# P2U Optimal SFO Pipeline

This folder keeps the realistic P2U optimization workflow modular:

1. `prepare_p2u_corridor_network.py`
   - Builds the transformer/source corridor candidate graph from the road graph.
   - Filters to the source-side road 2-edge component, prunes terminal corridors, and contracts degree-2 terminal chains for the ILP input.

2. `run_p2u_ilp_2edge.py`
   - Solves the minimum-length 2-edge-connected ILP backbone on the prepared candidate graph.
   - The ILP stage optimizes the backbone only.

3. `build_p2u_final_network.py`
   - Creates the final terminal-only network after the ILP.
   - Keeps graph nodes as transformers/sources only.
   - Keeps ILP backbone chains as edge load.
   - Adds radial tree attachments from non-backbone transformer terminals to the optimized backbone using road shortest-path distance.

4. `create_p2u_final_network_qgis.py`
   - Writes a temporary QGIS project for inspecting the final terminal graph.

5. `analyze_p2u_final_network.py`
   - Computes final topology, length, transformer/load coverage, and length-based failure-probability parameters.
   - Monte Carlo reliability is optional and off by default because the current terminal sampler is slow on this graph.

6. `compare_p2u_old_new_reliability.py`
   - Compares the original full P2U MV network with the optimized backbone-plus-trees network.
   - Reports topology and decomposed reliability terms: bridge/tree, non-bridge section edge-load, regular-chain internal risk, and optional generalized structural risk.
   - Generalized 3-edge structural chains are off by default because the current detector times out on the full P2U final graph.

The whole pipeline can be run from `run_p2u_optimal_pipeline.py`. During development, prefer running only the needed stages, for example:

```powershell
& C:\Users\rotem\anaconda3\envs\reliability\python.exe figures\optimal_sfo\build_p2u_final_network.py
& C:\Users\rotem\anaconda3\envs\reliability\python.exe figures\optimal_sfo\create_p2u_final_network_qgis.py
& C:\Users\rotem\anaconda3\envs\reliability\python.exe figures\optimal_sfo\analyze_p2u_final_network.py
& C:\Users\rotem\anaconda3\envs\reliability\python.exe figures\optimal_sfo\compare_p2u_old_new_reliability.py
```
