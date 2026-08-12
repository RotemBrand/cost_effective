from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.p2u_final_network import (  # noqa: E402
    FINAL_NETWORK_GPKG,
    FINAL_NETWORK_METADATA,
    OUTPUT_DIR,
    cycle_rank,
    graph_from_final_tables,
    read_final_network_tables,
)
from figures.p2u_full_mv_reliability import (  # noqa: E402
    aggregate_lv_consumers_to_transformers,
    build_full_mv_graph,
)
from indexes import GraphRel  # noqa: E402


OUTPUT_JSON = OUTPUT_DIR / "p2u_old_new_reliability_comparison.json"
OUTPUT_MD = OUTPUT_DIR / "p2u_old_new_reliability_comparison.md"


def _edge_key(u: Any, v: Any) -> tuple[Any, Any]:
    return tuple(sorted((u, v), key=repr))


def _prepare_old_graph() -> tuple[nx.Graph, list[Any], dict]:
    graph, mv_rows, switches = build_full_mv_graph()
    lv_assignment = aggregate_lv_consumers_to_transformers(graph)
    for _, _, data in graph.edges(data=True):
        data.setdefault("length", data.get("length_m", 1.0))
        data.setdefault("edge_weight", 0.0)
    return graph, list(graph.graph["sources"]), {
        "raw_mv_line_features": int(len(mv_rows)),
        "mv_switching_device_features": int(len(switches)),
        "lv_assignment": lv_assignment,
    }


def _prepare_new_graph(gpkg: Path, metadata_json: Path) -> tuple[nx.Graph, list[Any], dict]:
    metadata = json.loads(metadata_json.read_text(encoding="utf-8")) if metadata_json.exists() else {}
    if "output_gpkg" in metadata:
        gpkg = Path(metadata["output_gpkg"])
    nodes, edges, backbone_edges, tree_edges = read_final_network_tables(gpkg)
    graph = graph_from_final_tables(nodes, edges)
    sources = list(nodes[nodes["source_count"].astype(float) > 0]["terminal_id"].astype(str))
    return graph, sources, {
        "metadata": metadata,
        "backbone_edges": int(len(backbone_edges)),
        "tree_attachment_edges": int(len(tree_edges)),
        "tree_attachment_physical_road_union_length_m": float(
            metadata.get("tree_attachment_physical_road_union_length_m", 0.0)
        ),
    }


def _source_contracted_simple(graph: nx.Graph, sources: list[Any]) -> nx.Graph:
    source_set = set(sources)
    contracted = nx.Graph()
    for node, data in graph.nodes(data=True):
        target = "__SOURCE__" if node in source_set else node
        if target not in contracted:
            contracted.add_node(target, **data)
        else:
            contracted.nodes[target]["weight"] = float(contracted.nodes[target].get("weight", 0.0)) + float(
                data.get("weight", 0.0)
            )
    for u, v, data in graph.edges(data=True):
        cu = "__SOURCE__" if u in source_set else u
        cv = "__SOURCE__" if v in source_set else v
        if cu == cv:
            continue
        length = float(data.get("length", data.get("length_m", 1.0)))
        if contracted.has_edge(cu, cv) and contracted.edges[cu, cv].get("length", length) <= length:
            continue
        contracted.add_edge(cu, cv, **data)
    return contracted


def _graph_rel(graph: nx.Graph, sources: list[Any]) -> GraphRel:
    node_weights = {node: float(data.get("weight", 0.0)) for node, data in graph.nodes(data=True)}
    return GraphRel(graph, nodes_weight=node_weights, sources=sources)


def _poly_coeffs(poly) -> dict[str, float]:
    return {
        "p0": float(poly[0]),
        "p1": float(poly[1]),
        "p2": float(poly[2]),
    }


def _risk_terms_dict(terms, *, p_mean: float) -> dict:
    result = {
        "total_weight": float(terms.total_weight),
        "output": terms.output,
        "mean_edge_failure_prob": terms.mean_edge_failure_prob,
        "edge_failure_rate": terms.edge_failure_rate,
        "mean_actual_edge_probability": terms.mean_actual_edge_probability,
    }
    if terms.output == "float":
        result.update(
            {
                "tree": float(terms.tree),
                "nonbridge_section": float(terms.nonbridge_section),
                "internal_regular_chains": float(terms.internal),
                "structural_generalized_chains": float(terms.structural),
                "total": float(terms.total),
            }
        )
        return result

    parts = {
        "tree": terms.tree,
        "nonbridge_section": terms.nonbridge_section,
        "internal_regular_chains": terms.internal,
        "structural_generalized_chains": terms.structural,
        "total": terms.total,
    }
    for name, poly in parts.items():
        coeffs = _poly_coeffs(poly)
        coeffs["value_at_p_mean"] = float(poly(p_mean))
        result[name] = coeffs
    return result


def _topology_dict(
    *,
    graph: nx.Graph,
    sources: list[Any],
    graph_rel: GraphRel,
    decomp,
    extra: dict,
) -> dict:
    source_graph = decomp.source_graph
    source_contracted = _source_contracted_simple(graph, sources)
    total_node_weight = sum(float(data.get("weight", 0.0)) for _, data in graph.nodes(data=True))
    total_edge_weight = sum(float(data.get("edge_weight", 0.0)) for _, _, data in graph.edges(data=True))
    tie_count = sum(1 for _, _, data in graph.edges(data=True) if data.get("is_tie", False))
    switch_count = sum(1 for _, _, data in graph.edges(data=True) if data.get("is_switch", False))
    return {
        "raw_nodes": int(graph.number_of_nodes()),
        "raw_edges": int(graph.number_of_edges()),
        "raw_components": int(nx.number_connected_components(graph)),
        "sources": int(len(sources)),
        "raw_cycle_rank": int(cycle_rank(graph)),
        "tie_edges": int(tie_count),
        "switch_edges": int(switch_count),
        "total_length_m": float(sum(float(data.get("length", data.get("length_m", 0.0))) for _, _, data in graph.edges(data=True))),
        "node_weight": float(total_node_weight),
        "edge_weight": float(total_edge_weight),
        "raw_total_weight": float(total_node_weight + total_edge_weight),
        "reliability_total_weight_after_source_contraction": float(decomp.total_weight),
        "graphrel_nodes": int(graph_rel.graph.number_of_nodes()),
        "graphrel_edges": int(graph_rel.graph.number_of_edges()),
        "source_contracted_nodes": int(source_contracted.number_of_nodes()),
        "source_contracted_edges": int(source_contracted.number_of_edges()),
        "source_contracted_components": int(nx.number_connected_components(source_contracted)),
        "source_contracted_cycle_rank_R": int(cycle_rank(source_contracted)),
        "source_contracted_bridges": int(len(list(nx.bridges(source_contracted))) if nx.is_connected(source_contracted) else 0),
        "decomposition_bridges": int(len(decomp.bridges)),
        "two_edge_components": int(len(decomp.two_edge_components)),
        "bridge_tree_nodes": int(decomp.bridge_tree.number_of_nodes()),
        "bridge_tree_edges": int(decomp.bridge_tree.number_of_edges()),
        "structure_graph_nodes": int(decomp.structure_graph.number_of_nodes()),
        "structure_graph_edges": int(decomp.structure_graph.number_of_edges()),
        "regular_chains": int(len(decomp.regular_chains)),
        "regular_chain_total_length": float(sum(chain.length for chain in decomp.regular_chains)),
        "regular_chain_total_weight": float(sum(chain.total_weight for chain in decomp.regular_chains)),
        "three_edge_macro_nodes": int(decomp.three_edge_macro_graph.number_of_nodes()),
        "three_edge_macro_edges": int(decomp.three_edge_macro_graph.number_of_edges()),
        "generalized_chains": int(len(decomp.generalized_chains)),
        **extra,
    }


def _analyze_one(
    *,
    name: str,
    graph: nx.Graph,
    sources: list[Any],
    extra: dict,
    p_mean: float,
    generalized: bool,
    generalized_method: str,
) -> dict:
    t0 = time.perf_counter()
    graph_rel = _graph_rel(graph, sources)
    decomp = graph_rel.decompose(
        include_generalized_chains=generalized,
        generalized_component_method=generalized_method,
    )
    decompose_runtime = time.perf_counter() - t0

    float_terms = decomp.switch_risk_terms(
        mean_edge_failure_prob=p_mean,
        length_attr="length",
        output="float",
    )
    poly_terms = decomp.switch_risk_terms(
        mean_edge_failure_prob=p_mean,
        length_attr="length",
        output="poly",
    )
    return {
        "name": name,
        "generalized_structural_requested": bool(generalized),
        "generalized_component_method": generalized_method if generalized else None,
        "decomposition_runtime_seconds": float(decompose_runtime),
        "topology": _topology_dict(
            graph=graph,
            sources=sources,
            graph_rel=graph_rel,
            decomp=decomp,
            extra=extra,
        ),
        "risk_float": _risk_terms_dict(float_terms, p_mean=p_mean),
        "risk_poly": _risk_terms_dict(poly_terms, p_mean=p_mean),
    }


def compare_old_new(
    *,
    final_gpkg: Path = FINAL_NETWORK_GPKG,
    final_metadata: Path = FINAL_NETWORK_METADATA,
    p_mean: float = 5e-4,
    generalized: bool = False,
    generalized_method: str = "projection",
    output_json: Path = OUTPUT_JSON,
    output_md: Path = OUTPUT_MD,
) -> dict:
    old_graph, old_sources, old_extra = _prepare_old_graph()
    new_graph, new_sources, new_extra = _prepare_new_graph(final_gpkg, final_metadata)

    old = _analyze_one(
        name="original_full_p2u_mv",
        graph=old_graph,
        sources=old_sources,
        extra=old_extra,
        p_mean=p_mean,
        generalized=generalized,
        generalized_method=generalized_method,
    )
    new = _analyze_one(
        name="optimized_backbone_with_trees",
        graph=new_graph,
        sources=new_sources,
        extra=new_extra,
        p_mean=p_mean,
        generalized=generalized,
        generalized_method=generalized_method,
    )
    comparison = {
        "p_mean": float(p_mean),
        "generalized_structural_terms_computed": bool(generalized),
        "generalized_note": (
            "Generalized structural chains were computed with the selected 3-edge component method."
            if generalized
            else (
                "Generalized structural chains are off by default. With generalized=false, "
                "structural_generalized_chains is zero and the O(p^2) comparison includes "
                "tree p2, regular internal chains, and nonbridge section p1 only."
            )
        ),
        "old": old,
        "new": new,
        "delta_new_minus_old": _delta(old, new),
    }
    output_json.write_text(json.dumps(comparison, indent=2), encoding="utf-8")
    output_md.write_text(_markdown(comparison), encoding="utf-8")
    return comparison


def _get_path(data: dict, path: list[str], default=0.0):
    cur = data
    for key in path:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


def _delta(old: dict, new: dict) -> dict:
    paths = {
        "R": ["topology", "source_contracted_cycle_rank_R"],
        "bridges": ["topology", "source_contracted_bridges"],
        "length_m": ["topology", "total_length_m"],
        "total_weight": ["topology", "reliability_total_weight_after_source_contraction"],
        "float_total_risk": ["risk_float", "total"],
        "float_tree_risk": ["risk_float", "tree"],
        "float_nonbridge_section_risk": ["risk_float", "nonbridge_section"],
        "float_internal_regular_chain_risk": ["risk_float", "internal_regular_chains"],
        "float_structural_generalized_chain_risk": ["risk_float", "structural_generalized_chains"],
        "poly_total_p1": ["risk_poly", "total", "p1"],
        "poly_total_p2": ["risk_poly", "total", "p2"],
        "poly_total_value_at_p_mean": ["risk_poly", "total", "value_at_p_mean"],
    }
    delta = {}
    for name, path in paths.items():
        old_value = float(_get_path(old, path, 0.0))
        new_value = float(_get_path(new, path, 0.0))
        delta[name] = {
            "old": old_value,
            "new": new_value,
            "delta": new_value - old_value,
            "ratio_new_over_old": new_value / old_value if old_value != 0 else None,
        }
    return delta


def _fmt(value: Any, digits: int = 6) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return str(value)
    try:
        return f"{float(value):.{digits}g}"
    except (TypeError, ValueError):
        return str(value)


def _markdown(comparison: dict) -> str:
    old = comparison["old"]
    new = comparison["new"]
    rows = []
    metrics = [
        ("Source-contracted R", ["topology", "source_contracted_cycle_rank_R"]),
        ("Source-contracted bridges", ["topology", "source_contracted_bridges"]),
        ("Raw nodes", ["topology", "raw_nodes"]),
        ("Raw edges", ["topology", "raw_edges"]),
        ("Total length km", ["topology", "total_length_m"], 1000.0),
        ("Reliability weight kW", ["topology", "reliability_total_weight_after_source_contraction"]),
        ("2-edge components", ["topology", "two_edge_components"]),
        ("Structure graph nodes", ["topology", "structure_graph_nodes"]),
        ("Structure graph edges", ["topology", "structure_graph_edges"]),
        ("Regular chains", ["topology", "regular_chains"]),
        ("Float total risk", ["risk_float", "total"]),
        ("Float tree risk", ["risk_float", "tree"]),
        ("Float nonbridge section risk", ["risk_float", "nonbridge_section"]),
        ("Float internal regular-chain risk", ["risk_float", "internal_regular_chains"]),
        ("Float structural generalized-chain risk", ["risk_float", "structural_generalized_chains"]),
        ("Poly total p1 coeff", ["risk_poly", "total", "p1"]),
        ("Poly total p2 coeff", ["risk_poly", "total", "p2"]),
        ("Poly value at p_mean", ["risk_poly", "total", "value_at_p_mean"]),
    ]
    for metric in metrics:
        label, path = metric[0], metric[1]
        scale = metric[2] if len(metric) > 2 else 1.0
        old_value = float(_get_path(old, path, 0.0)) / scale
        new_value = float(_get_path(new, path, 0.0)) / scale
        rows.append(
            f"| {label} | `{_fmt(old_value)}` | `{_fmt(new_value)}` | `{_fmt(new_value - old_value)}` |"
        )

    return "\n".join(
        [
            "# P2U Old-New Reliability Comparison",
            "",
            f"- Mean edge failure probability target: `{comparison['p_mean']}`",
            f"- Generalized structural chains computed: `{comparison['generalized_structural_terms_computed']}`",
            "",
            comparison["generalized_note"],
            "",
            "| Metric | Original P2U MV | Optimized backbone + trees | New - old |",
            "|---|---:|---:|---:|",
            *rows,
            "",
            "## Decomposition Meaning",
            "",
            "- `tree`: bridge / 1-connected contribution after source contraction.",
            "- `nonbridge_section`: first-order risk from load stored on non-bridge section edges.",
            "- `internal_regular_chains`: second-order risk from two cuts inside regular chains.",
            "- `structural_generalized_chains`: second-order generalized-chain structural term. It is zero when `--generalized` is not used.",
            "",
            "## Runtime",
            "",
            f"- Original decomposition: `{old['decomposition_runtime_seconds']:.2f}` s",
            f"- New decomposition: `{new['decomposition_runtime_seconds']:.2f}` s",
        ]
    ) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare original and optimized P2U reliability decompositions.")
    parser.add_argument("--final-gpkg", type=Path, default=FINAL_NETWORK_GPKG)
    parser.add_argument("--final-metadata", type=Path, default=FINAL_NETWORK_METADATA)
    parser.add_argument("--p-mean", type=float, default=5e-4)
    parser.add_argument("--generalized", action="store_true", help="Also compute generalized 3-edge structural chains.")
    parser.add_argument("--generalized-method", choices=["projection", "networkx"], default="projection")
    parser.add_argument("--output-json", type=Path, default=OUTPUT_JSON)
    parser.add_argument("--output-md", type=Path, default=OUTPUT_MD)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    comparison = compare_old_new(
        final_gpkg=args.final_gpkg,
        final_metadata=args.final_metadata,
        p_mean=args.p_mean,
        generalized=args.generalized,
        generalized_method=args.generalized_method,
        output_json=args.output_json,
        output_md=args.output_md,
    )
    print(json.dumps(comparison, indent=2))


if __name__ == "__main__":
    main()
