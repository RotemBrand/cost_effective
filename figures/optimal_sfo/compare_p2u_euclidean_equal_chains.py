from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import geopandas as gpd
import networkx as nx

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from figures.optimal_sfo.compare_p2u_old_new_reliability import (  # noqa: E402
    _analyze_one,
    _delta,
    _fmt,
    _get_path,
    _prepare_old_graph,
)
from figures.optimal_sfo.run_p2u_euclidean_equal_chains import (  # noqa: E402
    OUTPUT_GPKG as DEFAULT_EUCLIDEAN_GPKG,
    OUTPUT_DIR as EUCLIDEAN_OUTPUT_DIR,
    SUMMARY_JSON as DEFAULT_EUCLIDEAN_SUMMARY,
)


OUTPUT_JSON = EUCLIDEAN_OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_reliability_comparison.json"
OUTPUT_MD = EUCLIDEAN_OUTPUT_DIR / "p2u_euclidean_equal_chains_Rmax50_reliability_comparison.md"


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def _euclidean_graph_from_gpkg(
    gpkg: Path,
    *,
    original_weight_by_terminal: dict[str, float],
) -> tuple[nx.Graph, list[str], dict[str, Any]]:
    node_layer, edge_layer = _detect_terminal_edge_layers(gpkg)
    nodes = gpd.read_file(gpkg, layer=node_layer)
    edges = gpd.read_file(gpkg, layer=edge_layer)

    graph = nx.Graph()
    missing_weights: list[str] = []
    for _, row in nodes.iterrows():
        terminal_id = str(row.terminal_id)
        weight = float(original_weight_by_terminal.get(terminal_id, 0.0))
        if str(row.kind) == "transformer" and terminal_id not in original_weight_by_terminal:
            missing_weights.append(terminal_id)
        graph.add_node(
            terminal_id,
            weight=weight,
            kind=str(row.kind),
            is_source=str(row.kind) == "source",
            source_count=int(row.source_count),
            transformer_count=int(row.transformer_count),
            size_kva=float(row.size_kva),
            nominal_voltage_kv=float(row.nominal_voltage_kv),
        )

    for _, row in edges.iterrows():
        u = str(row.terminal_a)
        v = str(row.terminal_b)
        length = float(row.length_m)
        graph.add_edge(
            u,
            v,
            length=length,
            length_m=length,
            edge_weight=0.0,
            is_tie=bool(row.is_tie),
            is_switch=bool(row.is_switch),
            normally_closed=bool(row.normally_closed),
        )

    sources = list(nodes.loc[nodes["kind"] == "source", "terminal_id"].astype(str))
    if missing_weights:
        raise ValueError(f"missing original demand weights for {len(missing_weights)} transformer terminals")
    if not sources:
        raise ValueError("Euclidean graph has no source terminals")
    return graph, sources, {
        "input_gpkg": str(gpkg),
        "node_layer": node_layer,
        "edge_layer": edge_layer,
        "euclidean_transformers": int((nodes["kind"] == "transformer").sum()),
        "euclidean_sources": int((nodes["kind"] == "source").sum()),
    }


def _detect_terminal_edge_layers(gpkg: Path) -> tuple[str, str]:
    import pyogrio

    layers = set(pyogrio.list_layers(gpkg)[:, 0])
    candidates = [
        ("euclidean_nodes", "euclidean_edges"),
        ("direct_ilp_nodes", "direct_ilp_edges"),
    ]
    for node_layer, edge_layer in candidates:
        if node_layer in layers and edge_layer in layers:
            return node_layer, edge_layer
    raise ValueError(f"no recognized terminal/edge layer pair in {gpkg}; found {sorted(layers)}")


def _comparison_indexes(original: dict[str, Any], euclidean: dict[str, Any]) -> dict[str, Any]:
    original_length = float(_get_path(original, ["topology", "total_length_m"], 0.0))
    euclidean_length = float(_get_path(euclidean, ["topology", "total_length_m"], 0.0))
    original_physical_r = float(_get_path(original, ["topology", "raw_cycle_rank"], 0.0))
    euclidean_physical_r = float(_get_path(euclidean, ["topology", "raw_cycle_rank"], 0.0))
    original_theory_r = float(_get_path(original, ["topology", "source_contracted_cycle_rank_R"], 0.0))
    euclidean_theory_r = float(_get_path(euclidean, ["topology", "source_contracted_cycle_rank_R"], 0.0))
    original_risk = float(_get_path(original, ["risk_float", "total"], 0.0))
    euclidean_risk = float(_get_path(euclidean, ["risk_float", "total"], 0.0))

    return {
        "Z_w_length_saving": 1.0 - euclidean_length / original_length if original_length else None,
        "Z_R_physical_cycle_rank_reduction": (
            1.0 - euclidean_physical_r / original_physical_r if original_physical_r else None
        ),
        "Z_R_source_contracted_reduction": (
            1.0 - euclidean_theory_r / original_theory_r if original_theory_r else None
        ),
        "Z_F_reliability_risk_reduction": 1.0 - euclidean_risk / original_risk if original_risk else None,
        "length_ratio_euclidean_over_original": euclidean_length / original_length if original_length else None,
        "risk_ratio_euclidean_over_original": euclidean_risk / original_risk if original_risk else None,
    }


def compare_euclidean_equal_chains(
    *,
    euclidean_gpkg: Path = DEFAULT_EUCLIDEAN_GPKG,
    euclidean_summary: Path = DEFAULT_EUCLIDEAN_SUMMARY,
    case_label: str = "euclidean_equal_chain",
    p_mean: float = 5e-4,
    generalized_method: str = "projection",
    output_json: Path = OUTPUT_JSON,
    output_md: Path = OUTPUT_MD,
) -> dict[str, Any]:
    old_graph, old_sources, old_extra = _prepare_old_graph()
    original_weight_by_terminal = {
        f"T:{node}": float(data.get("weight", 0.0))
        for node, data in old_graph.nodes(data=True)
    }
    euclidean_graph, euclidean_sources, euclidean_extra = _euclidean_graph_from_gpkg(
        euclidean_gpkg,
        original_weight_by_terminal=original_weight_by_terminal,
    )
    euclidean_extra["euclidean_summary"] = _read_json(euclidean_summary)

    t0 = time.perf_counter()
    original = _analyze_one(
        name="original_full_p2u_mv",
        graph=old_graph,
        sources=old_sources,
        extra=old_extra,
        p_mean=p_mean,
        generalized=True,
        generalized_method=generalized_method,
    )
    euclidean = _analyze_one(
        name=case_label,
        graph=euclidean_graph,
        sources=euclidean_sources,
        extra=euclidean_extra,
        p_mean=p_mean,
        generalized=True,
        generalized_method=generalized_method,
    )
    comparison = {
        "p_mean": float(p_mean),
        "failure_probability_convention": (
            "edge probabilities are length-scaled with mean_edge_failure_prob=p_mean, "
            "using length_attr='length'"
        ),
        "note": (
            "This compares the first-stage Euclidean equal-chain graph against the original P2U MV graph. "
            "The Euclidean graph uses straight-line edge lengths and original LV-assigned transformer demand weights; "
            "it is not yet embedded onto the street network."
        ),
        "runtime_seconds": float(time.perf_counter() - t0),
        "original": original,
        "euclidean": euclidean,
        "delta_euclidean_minus_original": _delta(original, euclidean),
        "indexes": _comparison_indexes(original, euclidean),
    }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(comparison, indent=2), encoding="utf-8")
    output_md.write_text(_markdown(comparison), encoding="utf-8")
    return comparison


def _markdown(comparison: dict[str, Any]) -> str:
    original = comparison["original"]
    euclidean = comparison["euclidean"]
    indexes = comparison["indexes"]

    rows = [
        ("Total cost / length km", ["topology", "total_length_m"], 1000.0),
        ("Physical cycle rank R", ["topology", "raw_cycle_rank"], 1.0),
        ("Source-contracted R", ["topology", "source_contracted_cycle_rank_R"], 1.0),
        ("Tie / normally-open edges", ["topology", "tie_edges"], 1.0),
        ("Reliability weight kW", ["topology", "reliability_total_weight_after_source_contraction"], 1.0),
        ("Total risk F", ["risk_float", "total"], 1.0),
        ("Tree risk", ["risk_float", "tree"], 1.0),
        ("Nonbridge section risk", ["risk_float", "nonbridge_section"], 1.0),
        ("Internal regular-chain risk", ["risk_float", "internal_regular_chains"], 1.0),
        ("Generalized structural risk", ["risk_float", "structural_generalized_chains"], 1.0),
        ("Total p1 coefficient", ["risk_poly", "total", "p1"], 1.0),
        ("Total p2 coefficient", ["risk_poly", "total", "p2"], 1.0),
    ]
    table = []
    for label, path, scale in rows:
        original_value = float(_get_path(original, path, 0.0)) / scale
        euclidean_value = float(_get_path(euclidean, path, 0.0)) / scale
        table.append(
            f"| {label} | `{_fmt(original_value)}` | `{_fmt(euclidean_value)}` | `{_fmt(euclidean_value - original_value)}` |"
        )

    return "\n".join(
        [
            "# P2U Euclidean Equal-Chain Reliability Comparison",
            "",
            f"- Mean edge failure probability target: `{comparison['p_mean']}`",
            f"- Failure convention: {comparison['failure_probability_convention']}",
            f"- Total runtime: `{comparison['runtime_seconds']:.2f}` s",
            "",
            comparison["note"],
            "",
            "## Main Indexes",
            "",
            f"- `Z_w = 1 - W_new / W_original`: `{_fmt(indexes['Z_w_length_saving'])}`",
            f"- `Z_R`, physical cycle rank: `{_fmt(indexes['Z_R_physical_cycle_rank_reduction'])}`",
            f"- `Z_R`, source-contracted theory graph: `{_fmt(indexes['Z_R_source_contracted_reduction'])}`",
            f"- `Z_F = 1 - F_new / F_original`: `{_fmt(indexes['Z_F_reliability_risk_reduction'])}`",
            "",
            f"| Metric | Original P2U MV | {euclidean['name']} | New - original |",
            "|---|---:|---:|---:|",
            *table,
            "",
            "## Interpretation Notes",
            "",
            (
                "- The Euclidean construction reports raw physical `R`; reliability contracts all 15 source nodes "
                "into one source, so the source-contracted theory `R` can be larger."
            ),
            "- The cost here is straight-line Euclidean length. Road embedding will change the cost and should be treated as the next-stage feasibility/cost test.",
            "- The Euclidean equal-chain graph has no bridge/tree risk in this calculation; the remaining risk is the exact `p^2` regular-chain internal term.",
            "",
            "## Decomposition Runtime",
            "",
            f"- Original P2U MV: `{original['decomposition_runtime_seconds']:.2f}` s",
            f"- {euclidean['name']}: `{euclidean['decomposition_runtime_seconds']:.2f}` s",
        ]
    ) + "\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare P2U original MV reliability with Euclidean equal-chain R50.")
    parser.add_argument("--euclidean-gpkg", type=Path, default=DEFAULT_EUCLIDEAN_GPKG)
    parser.add_argument("--euclidean-summary", type=Path, default=DEFAULT_EUCLIDEAN_SUMMARY)
    parser.add_argument("--case-label", default="euclidean_equal_chain")
    parser.add_argument("--p-mean", type=float, default=5e-4)
    parser.add_argument("--generalized-method", choices=["projection", "networkx"], default="projection")
    parser.add_argument("--output-json", type=Path, default=OUTPUT_JSON)
    parser.add_argument("--output-md", type=Path, default=OUTPUT_MD)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    comparison = compare_euclidean_equal_chains(
        euclidean_gpkg=args.euclidean_gpkg,
        euclidean_summary=args.euclidean_summary,
        case_label=args.case_label,
        p_mean=args.p_mean,
        generalized_method=args.generalized_method,
        output_json=args.output_json,
        output_md=args.output_md,
    )
    print(json.dumps(comparison, indent=2))


if __name__ == "__main__":
    main()
