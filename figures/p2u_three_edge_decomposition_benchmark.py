from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from indexes import GraphRel
from figures.p2u_full_mv_reliability import (
    OUT_DIR,
    aggregate_lv_consumers_to_transformers,
    build_full_mv_graph,
)


def run_benchmark(*, method: str, include_generalized: bool) -> dict:
    graph, _, _ = build_full_mv_graph()
    lv_assignment = aggregate_lv_consumers_to_transformers(graph)
    sources = graph.graph["sources"]
    weights = {node: float(data.get("weight", 0.0)) for node, data in graph.nodes(data=True)}

    t0 = time.perf_counter()
    graph_rel = GraphRel(graph, nodes_weight=weights, sources=sources)
    graphrel_seconds = time.perf_counter() - t0

    t1 = time.perf_counter()
    decomp = graph_rel.decompose(
        include_generalized_chains=include_generalized,
        generalized_component_method=method,
    )
    decompose_seconds = time.perf_counter() - t1

    generalized_lengths = sorted((float(chain.length) for chain in decomp.generalized_chains), reverse=True)
    regular_lengths = sorted((float(chain.length) for chain in decomp.regular_chains), reverse=True)
    return {
        "dataset": "SMART-DS SFO P2U full MV network",
        "method": method,
        "include_generalized_chains": include_generalized,
        "graph": {
            "nodes": graph.number_of_nodes(),
            "edges": graph.number_of_edges(),
            "sources": len(sources),
            "total_weight_kw": sum(weights.values()),
        },
        "graphrel": {
            "nodes": graph_rel.graph.number_of_nodes(),
            "edges": graph_rel.graph.number_of_edges(),
            "init_seconds": graphrel_seconds,
        },
        "decomposition": {
            "seconds": decompose_seconds,
            "bridges": len(decomp.bridges),
            "two_edge_components": len(decomp.two_edge_components),
            "structure_graph_nodes": decomp.structure_graph.number_of_nodes(),
            "structure_graph_edges": decomp.structure_graph.number_of_edges(),
            "regular_chains": len(decomp.regular_chains),
            "three_edge_macro_graph_nodes": decomp.three_edge_macro_graph.number_of_nodes(),
            "three_edge_macro_graph_edges": decomp.three_edge_macro_graph.number_of_edges(),
            "generalized_chains": len(decomp.generalized_chains),
            "parallel_macro_edges": sum(
                1
                for _, _, data in decomp.three_edge_macro_graph.edges(data=True)
                if int(data.get("parallel_macro_edge_count", 1)) > 1
            ),
            "largest_regular_chain_lengths_m": regular_lengths[:10],
            "largest_generalized_chain_lengths_m": generalized_lengths[:10],
        },
        "lv_assignment": lv_assignment,
    }


def write_outputs(metrics: dict, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    suffix = metrics["method"] if metrics["include_generalized_chains"] else f"{metrics['method']}_no_generalized"
    stem = f"p2u_three_edge_decomposition_{suffix}"
    json_path = output_dir / f"{stem}.json"
    md_path = output_dir / f"{stem}.md"
    json_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    graph = metrics["graph"]
    graphrel = metrics["graphrel"]
    decomp = metrics["decomposition"]
    md_path.write_text(
        "\n".join(
            [
                "# P2U 3-Edge Decomposition Benchmark",
                "",
                f"- Method: `{metrics['method']}`",
                f"- Include generalized chains: `{metrics['include_generalized_chains']}`",
                f"- Original graph: `{graph['nodes']}` nodes, `{graph['edges']}` edges, `{graph['sources']}` sources",
                f"- GraphRel graph: `{graphrel['nodes']}` nodes, `{graphrel['edges']}` edges",
                f"- GraphRel init runtime: `{graphrel['init_seconds']:.3f}` s",
                f"- Decomposition runtime: `{decomp['seconds']:.3f}` s",
                f"- Bridges: `{decomp['bridges']}`",
                f"- 2-edge components: `{decomp['two_edge_components']}`",
                f"- Structure graph: `{decomp['structure_graph_nodes']}` nodes, `{decomp['structure_graph_edges']}` edges",
                f"- Regular chains: `{decomp['regular_chains']}`",
                f"- 3-edge macro graph: `{decomp['three_edge_macro_graph_nodes']}` nodes, `{decomp['three_edge_macro_graph_edges']}` edges",
                f"- Generalized chains: `{decomp['generalized_chains']}`",
                f"- Aggregated parallel macro edges: `{decomp['parallel_macro_edges']}`",
                f"- Largest generalized chain lengths, m: `{decomp['largest_generalized_chain_lengths_m']}`",
                "",
                "Run command:",
                "",
                "```powershell",
                f"conda run -n reliability python figures\\p2u_three_edge_decomposition_benchmark.py --method {metrics['method']}",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark P2U 3-edge decomposition backends.")
    parser.add_argument("--method", choices=["networkx", "projection"], default="projection")
    parser.add_argument("--no-generalized", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    metrics = run_benchmark(method=args.method, include_generalized=not args.no_generalized)
    write_outputs(metrics, args.output_dir)
    print(json.dumps(metrics["decomposition"], indent=2))


if __name__ == "__main__":
    main()
