from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import geopandas as gpd
import gurobipy as gp
import networkx as nx
import pandas as pd
from gurobipy import GRB


ROOT = Path(__file__).resolve().parents[2]
INPUT_GPKG = ROOT / "outputs" / "optimal_sfo" / "p2u_terminal_corridors_road2_k10_3857.gpkg"
OUTPUT_DIR = ROOT / "outputs" / "optimal_sfo"
OUTPUT_GPKG = OUTPUT_DIR / "p2u_ilp_2edge_solution_3857.gpkg"
SUMMARY_JSON = OUTPUT_DIR / "p2u_ilp_2edge_solution_summary.json"
SUMMARY_MD = OUTPUT_DIR / "p2u_ilp_2edge_solution_summary.md"
SUPER_SOURCE = "__SOURCE__"


def _edge_key(u: str, v: str) -> tuple[str, str]:
    return tuple(sorted((str(u), str(v))))


def _terminal_road_node(terminal: str) -> str:
    return str(terminal).split(":", 1)[1] if ":" in str(terminal) else str(terminal)


def _raw_terminal_edge_lengths(input_gpkg: Path) -> dict[tuple[str, str], float]:
    try:
        raw_edges = gpd.read_file(input_gpkg, layer="raw_terminal_edges")
    except Exception:
        return {}
    lengths: dict[tuple[str, str], float] = {}
    for _, row in raw_edges.iterrows():
        key = _edge_key(str(row.terminal_a), str(row.terminal_b))
        length = float(row.length_m)
        if key not in lengths or length < lengths[key]:
            lengths[key] = length
    return lengths


def _chain_cost_metrics(row, raw_lengths: dict[tuple[str, str], float]) -> dict[str, float | str]:
    length_m = float(row.length_m)
    sequence = [token for token in str(row.terminal_sequence).split("|") if token]
    segment_lengths: list[float] = []
    for a, b in zip(sequence, sequence[1:]):
        segment = raw_lengths.get(_edge_key(a, b))
        if segment is not None:
            segment_lengths.append(float(segment))
    if not segment_lengths:
        n_segments = max(1, len(sequence) - 1)
        segment_lengths = [length_m / n_segments for _ in range(n_segments)]

    max_segment = max(segment_lengths) if segment_lengths else length_m
    has_internal_terminals = int(row.internal_terminal_count) > 0
    tree_service = max(0.0, length_m - max_segment) if has_internal_terminals else 0.0
    extra_backbone = length_m - tree_service
    return {
        "chain_total_length_m": length_m,
        "chain_max_segment_length_m": float(max_segment),
        "chain_tree_service_length_m": float(tree_service),
        "chain_extra_backbone_cost_m": float(extra_backbone),
        "terminal_sequence": str(row.terminal_sequence),
    }


def load_ilp_graph(input_gpkg: Path = INPUT_GPKG) -> tuple[nx.Graph, gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame, dict[str, set[tuple[str, str]]]]:
    edges = gpd.read_file(input_gpkg, layer="ilp_edges")
    raw_lengths = _raw_terminal_edge_lengths(input_gpkg)
    transformer_nodes = gpd.read_file(input_gpkg, layer="ilp_transformer_nodes")
    source_nodes = gpd.read_file(input_gpkg, layer="ilp_source_nodes")
    nodes = pd.concat([transformer_nodes, source_nodes], ignore_index=True)
    nodes = nodes.drop_duplicates(subset=["terminal_id"], keep="first")

    graph = nx.Graph()
    original_source_terminals = set(source_nodes["terminal_id"].astype(str))
    original_source_incident: dict[str, set[tuple[str, str]]] = {
        source: set() for source in original_source_terminals
    }
    for _, row in nodes.iterrows():
        terminal = str(row.terminal_id)
        node = SUPER_SOURCE if str(row.kind) in {"source", "source_transformer"} else terminal
        if node in graph:
            graph.nodes[node]["transformer_count"] += int(row.transformer_count)
            graph.nodes[node]["source_count"] += int(row.source_count)
            graph.nodes[node]["size_kva"] += float(row.size_kva)
            continue
        graph.add_node(
            node,
            kind="source" if node == SUPER_SOURCE else str(row.kind),
            road_node=SUPER_SOURCE if node == SUPER_SOURCE else str(row.road_node),
            transformer_count=int(row.transformer_count),
            source_count=int(row.source_count),
            size_kva=float(row.size_kva),
        )

    for edge_id, row in edges.reset_index(drop=True).iterrows():
        u_raw = str(row.terminal_a)
        v_raw = str(row.terminal_b)
        kind_u = str(row.terminal_kind_a)
        kind_v = str(row.terminal_kind_b)
        u = SUPER_SOURCE if kind_u in {"source", "source_transformer"} else u_raw
        v = SUPER_SOURCE if kind_v in {"source", "source_transformer"} else v_raw
        if u == v:
            continue
        length_m = float(row.length_m)
        if graph.has_edge(u, v) and graph.edges[u, v]["length_m"] <= length_m:
            continue
        graph.add_edge(
            u,
            v,
            edge_id=int(edge_id),
            original_terminal_a=u_raw,
            original_terminal_b=v_raw,
            length_m=length_m,
            edge_transformer_count=int(row.edge_transformer_count),
            edge_size_kva=float(row.edge_size_kva),
            **_chain_cost_metrics(row, raw_lengths),
        )
        if u_raw in original_source_incident:
            original_source_incident[u_raw].add(_edge_key(u, v))
        if v_raw in original_source_incident:
            original_source_incident[v_raw].add(_edge_key(u, v))
    return graph, edges, transformer_nodes, source_nodes, original_source_incident


def solve_min_2edge(
    graph: nx.Graph,
    *,
    original_source_incident: dict[str, set[tuple[str, str]]] | None,
    redundancy: int | None,
    max_redundancy: int | None,
    time_limit: float,
    mip_gap: float,
    threads: int,
    max_cut_rounds: int,
    cut_mode: str = "iterative",
    warm_start_edges: set[tuple[str, str]] | None = None,
    objective_mode: str = "min_length",
    coverage_attr: str = "edge_size_kva",
    coverage_tolerance: float = 1e-6,
    source_min_incident_edges: int = 1,
    chain_extra_cost_budget_m: float | None = None,
) -> tuple[nx.Graph | None, dict]:
    if redundancy is not None and max_redundancy is not None:
        raise ValueError("Use either exact redundancy or max_redundancy, not both")
    if cut_mode not in {"iterative", "callback"}:
        raise ValueError("cut_mode must be 'iterative' or 'callback'")
    if objective_mode not in {"min_length", "max_chain_demand_then_min_length", "max_chain_demand_under_cost"}:
        raise ValueError(
            "objective_mode must be 'min_length', 'max_chain_demand_then_min_length', "
            "or 'max_chain_demand_under_cost'"
        )
    if source_min_incident_edges < 0:
        raise ValueError("source_min_incident_edges must be nonnegative")

    nodes = sorted(graph.nodes)
    edges = sorted(_edge_key(u, v) for u, v in graph.edges)
    incident = {node: [] for node in nodes}
    for u, v in edges:
        incident[u].append((u, v))
        incident[v].append((u, v))

    model = gp.Model("p2u_min_2edge")
    model.Params.OutputFlag = 0
    model.Params.MIPFocus = 1
    model.Params.MIPGap = mip_gap
    model.Params.Threads = threads
    if time_limit > 0:
        model.Params.TimeLimit = time_limit
    if cut_mode == "callback":
        model.Params.LazyConstraints = 1

    x = model.addVars(edges, vtype=GRB.BINARY, name="x")
    warm_start_edges = {_edge_key(u, v) for u, v in (warm_start_edges or set())}
    warm_start_edges = {edge for edge in warm_start_edges if edge in x}
    if warm_start_edges:
        for edge in edges:
            x[edge].Start = 1.0 if edge in warm_start_edges else 0.0
    for node in nodes:
        model.addConstr(gp.quicksum(x[e] for e in incident[node]) >= 2, name=f"deg2[{node}]")

    source_connectivity_constraints = 0
    insufficient_source_degree_candidates: dict[str, int] = {}
    if original_source_incident:
        for source, source_edges in sorted(original_source_incident.items()):
            available_edges = [edge for edge in source_edges if edge in x]
            if not available_edges:
                continue
            enforced_min = min(int(source_min_incident_edges), len(available_edges))
            if len(available_edges) < int(source_min_incident_edges):
                insufficient_source_degree_candidates[source] = len(available_edges)
            model.addConstr(
                gp.quicksum(x[e] for e in available_edges) >= enforced_min,
                name=f"source_incident[{source}]",
            )
            source_connectivity_constraints += 1

    if redundancy is not None:
        model.addConstr(
            gp.quicksum(x[e] for e in edges) == len(nodes) - 1 + int(redundancy),
            name="edge_budget",
        )
    if max_redundancy is not None:
        model.addConstr(
            gp.quicksum(x[e] for e in edges) <= len(nodes) - 1 + int(max_redundancy),
            name="max_edge_budget",
        )

    length_expr = gp.quicksum(float(graph.edges[e]["length_m"]) * x[e] for e in edges)
    coverage_weights = {e: float(graph.edges[e].get(coverage_attr, 0.0) or 0.0) for e in edges}
    coverage_expr = gp.quicksum(coverage_weights[e] * x[e] for e in edges)
    chain_extra_costs = {e: float(graph.edges[e].get("chain_extra_backbone_cost_m", graph.edges[e]["length_m"]) or 0.0) for e in edges}
    chain_extra_cost_expr = gp.quicksum(chain_extra_costs[e] * x[e] for e in edges)
    if chain_extra_cost_budget_m is not None:
        model.addConstr(
            chain_extra_cost_expr <= float(chain_extra_cost_budget_m),
            name="chain_extra_cost_budget",
        )
    node_coverage_weight = sum(float(graph.nodes[node].get("size_kva", 0.0) or 0.0) for node in nodes)
    total_edge_coverage_weight = sum(coverage_weights.values())
    total_chain_tree_service_length_m = sum(
        float(graph.edges[e].get("chain_tree_service_length_m", 0.0) or 0.0) for e in edges
    )
    total_chain_extra_cost_m = sum(chain_extra_costs.values())
    objective_phases: list[dict] = []

    start = time.time()
    cut_round = 0
    total_added_cuts = 0
    selected_edges: list[tuple[str, str]] = []
    stop_reason = "max_cut_rounds"

    def violated_cutsets(selected: list[tuple[str, str]]) -> list[list[tuple[str, str]]]:
        solution = nx.Graph()
        solution.add_nodes_from(nodes)
        solution.add_edges_from(selected)

        cutsets: list[list[tuple[str, str]]] = []
        components = list(nx.connected_components(solution))
        if len(components) > 1:
            for comp in components:
                crossing = [e for e in edges if (e[0] in comp) ^ (e[1] in comp)]
                if crossing:
                    cutsets.append(crossing)
            return cutsets

        bridges = list(nx.bridges(solution))
        if not bridges:
            return cutsets
        bridge_keys = {_edge_key(u, v) for u, v in bridges}
        bridge_removed = solution.copy()
        bridge_removed.remove_edges_from(bridge_keys)
        for comp in nx.connected_components(bridge_removed):
            crossing = [e for e in edges if (e[0] in comp) ^ (e[1] in comp)]
            if crossing:
                cutsets.append(crossing)
        return cutsets

    def run_optimization_phase(phase_name: str) -> None:
        nonlocal cut_round, total_added_cuts, selected_edges, stop_reason
        phase_start = time.time()
        cuts_before = total_added_cuts
        rounds_before = cut_round
        if cut_mode == "callback":
            callback_stats = {"rounds": 0, "cuts": 0}

            def add_lazy_cuts(cb_model, where):
                if where != GRB.Callback.MIPSOL:
                    return
                vals = cb_model.cbGetSolution(x)
                incumbent_edges = [e for e in edges if vals[e] > 0.5]
                cutsets = violated_cutsets(incumbent_edges)
                if not cutsets:
                    return
                callback_stats["rounds"] += 1
                callback_stats["cuts"] += len(cutsets)
                for crossing in cutsets:
                    cb_model.cbLazy(gp.quicksum(x[e] for e in crossing) >= 2)

            if time_limit > 0:
                remaining_time = time_limit - (time.time() - start)
                if remaining_time <= 0:
                    stop_reason = "global_time_limit"
                    objective_phases.append(
                        {
                            "phase": phase_name,
                            "status_name": "SKIPPED_TIME_LIMIT",
                            "runtime_s": 0.0,
                            "cut_rounds": 0,
                            "added_cuts": 0,
                        }
                    )
                    return
                model.Params.TimeLimit = remaining_time

            model.optimize(add_lazy_cuts)
            cut_round += int(callback_stats["rounds"])
            total_added_cuts += int(callback_stats["cuts"])
            if model.SolCount == 0:
                stop_reason = "no_solution"
            else:
                selected_edges = [e for e in edges if x[e].X > 0.5]
                if violated_cutsets(selected_edges):
                    stop_reason = "time_limit_with_violated_cuts" if model.Status == GRB.TIME_LIMIT else "violated_cuts_remaining"
                else:
                    stop_reason = "2edge_connected"
        else:
            while cut_round < max_cut_rounds:
                if time_limit > 0:
                    remaining_time = time_limit - (time.time() - start)
                    if remaining_time <= 0:
                        stop_reason = "global_time_limit"
                        break
                    model.Params.TimeLimit = remaining_time
                cut_round += 1
                model.optimize()
                if model.SolCount == 0:
                    stop_reason = "no_solution"
                    break

                selected_edges = [e for e in edges if x[e].X > 0.5]
                cutsets = violated_cutsets(selected_edges)
                if not cutsets:
                    stop_reason = "2edge_connected"
                    break

                for crossing in cutsets:
                    model.addConstr(gp.quicksum(x[e] for e in crossing) >= 2)

                total_added_cuts += len(cutsets)
                model.update()

        objective_phases.append(
            {
                "phase": phase_name,
                "status": int(model.Status),
                "status_name": _gurobi_status_name(model.Status),
                "stop_reason": stop_reason,
                "runtime_s": time.time() - phase_start,
                "cut_rounds": cut_round - rounds_before,
                "added_cuts": total_added_cuts - cuts_before,
                "solution_count": int(model.SolCount),
                "objective_value": float(model.ObjVal) if model.SolCount else None,
                "best_bound": float(model.ObjBound) if model.SolCount else None,
                "mip_gap": float(model.MIPGap) if model.SolCount else None,
                "selected_edge_coverage_weight": float(sum(coverage_weights[e] for e in selected_edges)),
                "selected_length_m": float(sum(float(graph.edges[e]["length_m"]) for e in selected_edges)),
                "selected_chain_extra_cost_m": float(sum(chain_extra_costs[e] for e in selected_edges)),
            }
        )

    if objective_mode in {"max_chain_demand_then_min_length", "max_chain_demand_under_cost"}:
        model.setObjective(coverage_expr, GRB.MAXIMIZE)
        run_optimization_phase("maximize_chain_demand")
        if model.SolCount > 0:
            selected_coverage = float(sum(coverage_weights[e] for e in selected_edges))
            remaining_time = time_limit - (time.time() - start) if time_limit > 0 else None
            if remaining_time is None or remaining_time > 1.0:
                model.addConstr(
                    coverage_expr >= max(0.0, selected_coverage - float(coverage_tolerance)),
                    name="selected_chain_demand_floor",
                )
                model.setObjective(length_expr, GRB.MINIMIZE)
                model.update()
                run_optimization_phase("minimize_length_at_chain_demand")
            else:
                objective_phases.append(
                    {
                        "phase": "minimize_length_at_chain_demand",
                        "status_name": "SKIPPED_TIME_LIMIT",
                        "runtime_s": 0.0,
                        "cut_rounds": 0,
                        "added_cuts": 0,
                        "selected_edge_coverage_weight": selected_coverage,
                        "selected_length_m": float(sum(float(graph.edges[e]["length_m"]) for e in selected_edges)),
                        "selected_chain_extra_cost_m": float(sum(chain_extra_costs[e] for e in selected_edges)),
                    }
                )
    else:
        model.setObjective(length_expr, GRB.MINIMIZE)
        run_optimization_phase("minimize_length")

    elapsed = time.time() - start
    selected_length_m = float(sum(float(graph.edges[e]["length_m"]) for e in selected_edges))
    selected_edge_coverage_weight = float(sum(coverage_weights[e] for e in selected_edges))
    selected_chain_extra_cost_m = float(sum(chain_extra_costs[e] for e in selected_edges))
    selected_total_coverage_weight = float(node_coverage_weight + selected_edge_coverage_weight)
    if model.SolCount == 0:
        return None, {
            "status": int(model.Status),
            "status_name": _gurobi_status_name(model.Status),
            "stop_reason": stop_reason,
            "runtime_s": elapsed,
            "cut_rounds": cut_round,
            "added_cuts": total_added_cuts,
            "input_nodes": graph.number_of_nodes(),
            "input_edges": graph.number_of_edges(),
            "redundancy_constraint": redundancy,
            "max_redundancy_constraint": max_redundancy,
            "time_limit_s": time_limit,
            "threads": threads,
            "requested_mip_gap": mip_gap,
            "cut_mode": cut_mode,
            "source_connectivity_constraints": source_connectivity_constraints,
            "source_min_incident_edges": int(source_min_incident_edges),
            "insufficient_source_degree_candidates": insufficient_source_degree_candidates,
            "warm_start_edges": len(warm_start_edges),
            "objective_mode": objective_mode,
            "coverage_attr": coverage_attr,
            "chain_extra_cost_budget_m": None if chain_extra_cost_budget_m is None else float(chain_extra_cost_budget_m),
            "objective_phases": objective_phases,
            "best_bound_phase": objective_phases[-1]["phase"] if objective_phases else None,
            "node_coverage_weight": float(node_coverage_weight),
            "total_edge_coverage_weight": float(total_edge_coverage_weight),
            "total_chain_tree_service_length_m": float(total_chain_tree_service_length_m),
            "total_chain_extra_cost_m": float(total_chain_extra_cost_m),
        }

    solution = nx.Graph()
    solution.add_nodes_from((node, graph.nodes[node]) for node in nodes)
    for u, v in selected_edges:
        solution.add_edge(u, v, **graph.edges[u, v])

    bridges = list(nx.bridges(solution)) if nx.is_connected(solution) else []
    summary = {
        "status": int(model.Status),
        "status_name": _gurobi_status_name(model.Status),
        "stop_reason": stop_reason,
        "runtime_s": elapsed,
        "cut_rounds": cut_round,
        "added_cuts": total_added_cuts,
        "objective_length_m": selected_length_m,
        "best_bound": float(model.ObjBound),
        "mip_gap": float(model.MIPGap) if model.SolCount else None,
        "input_nodes": graph.number_of_nodes(),
        "input_edges": graph.number_of_edges(),
        "solution_nodes": solution.number_of_nodes(),
        "solution_edges": solution.number_of_edges(),
        "solution_cycle_rank": solution.number_of_edges() - solution.number_of_nodes() + nx.number_connected_components(solution),
        "solution_connected": nx.is_connected(solution),
        "solution_bridge_count": len(bridges),
        "solution_is_2edge_connected": nx.is_connected(solution) and len(bridges) == 0,
        "redundancy_constraint": redundancy,
        "max_redundancy_constraint": max_redundancy,
        "time_limit_s": time_limit,
        "threads": threads,
        "requested_mip_gap": mip_gap,
        "cut_mode": cut_mode,
        "source_connectivity_constraints": source_connectivity_constraints,
        "source_min_incident_edges": int(source_min_incident_edges),
        "insufficient_source_degree_candidates": insufficient_source_degree_candidates,
        "warm_start_edges": len(warm_start_edges),
        "objective_mode": objective_mode,
        "coverage_attr": coverage_attr,
        "chain_extra_cost_budget_m": None if chain_extra_cost_budget_m is None else float(chain_extra_cost_budget_m),
        "objective_phases": objective_phases,
        "best_bound_phase": objective_phases[-1]["phase"] if objective_phases else None,
        "node_coverage_weight": float(node_coverage_weight),
        "total_edge_coverage_weight": float(total_edge_coverage_weight),
        "total_chain_tree_service_length_m": float(total_chain_tree_service_length_m),
        "total_chain_extra_cost_m": float(total_chain_extra_cost_m),
        "selected_chain_extra_cost_m": selected_chain_extra_cost_m,
        "selected_edge_coverage_weight": selected_edge_coverage_weight,
        "selected_total_coverage_weight": selected_total_coverage_weight,
        "selected_edge_coverage_fraction": (
            selected_edge_coverage_weight / total_edge_coverage_weight if total_edge_coverage_weight > 0 else None
        ),
        "selected_total_coverage_fraction": (
            selected_total_coverage_weight / (node_coverage_weight + total_edge_coverage_weight)
            if node_coverage_weight + total_edge_coverage_weight > 0
            else None
        ),
        "selected_original_sources_with_incident_edge": _count_selected_sources(solution, original_source_incident),
        "original_sources_with_candidate_edges": len([s for s, es in (original_source_incident or {}).items() if es]),
    }
    return solution, summary


def _count_selected_sources(
    solution: nx.Graph,
    original_source_incident: dict[str, set[tuple[str, str]]] | None,
) -> int:
    if not original_source_incident:
        return 0
    selected = {_edge_key(u, v) for u, v in solution.edges}
    return sum(1 for edges in original_source_incident.values() if selected & edges)


def _gurobi_status_name(status: int) -> str:
    names = {
        GRB.LOADED: "LOADED",
        GRB.OPTIMAL: "OPTIMAL",
        GRB.INFEASIBLE: "INFEASIBLE",
        GRB.TIME_LIMIT: "TIME_LIMIT",
        GRB.INTERRUPTED: "INTERRUPTED",
        GRB.SUBOPTIMAL: "SUBOPTIMAL",
        GRB.INF_OR_UNBD: "INF_OR_UNBD",
        GRB.UNBOUNDED: "UNBOUNDED",
    }
    return names.get(status, f"STATUS_{status}")


def write_solution(
    solution: nx.Graph | None,
    summary: dict,
    input_edges: gpd.GeoDataFrame,
    transformer_nodes: gpd.GeoDataFrame,
    source_nodes: gpd.GeoDataFrame,
    *,
    output_gpkg: Path = OUTPUT_GPKG,
    summary_json: Path = SUMMARY_JSON,
    summary_md: Path = SUMMARY_MD,
) -> None:
    output_gpkg.parent.mkdir(parents=True, exist_ok=True)
    actual_output_gpkg = output_gpkg
    if output_gpkg.exists():
        try:
            output_gpkg.unlink()
        except PermissionError:
            actual_output_gpkg = output_gpkg.with_name(f"{output_gpkg.stem}_{int(time.time())}{output_gpkg.suffix}")
    summary["output_gpkg"] = str(actual_output_gpkg)

    if solution is not None:
        selected_ids = {int(data["edge_id"]) for _, _, data in solution.edges(data=True)}
        rows = input_edges.reset_index(drop=True).loc[lambda df: df.index.isin(selected_ids)].copy()
        rows.to_file(actual_output_gpkg, layer="solution_edges", driver="GPKG")

        selected_nodes = set(solution.nodes)
        transformer_nodes[transformer_nodes["terminal_id"].isin(selected_nodes)].to_file(
            actual_output_gpkg, layer="solution_transformer_nodes", driver="GPKG"
        )
        source_selection = source_nodes if SUPER_SOURCE in selected_nodes else source_nodes.iloc[0:0]
        source_selection.to_file(
            actual_output_gpkg, layer="solution_source_nodes", driver="GPKG"
        )

    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    summary_md.write_text(_summary_markdown(summary), encoding="utf-8")


def _summary_markdown(summary: dict) -> str:
    best_bound = summary.get("best_bound")
    best_bound_text = "" if best_bound is None else f"{float(best_bound):.3f}"
    lines = [
        "# P2U ILP 2-Edge Solution",
        "",
        f"- Status: `{summary.get('status_name')}`",
        f"- Stop reason: `{summary.get('stop_reason')}`",
        f"- Runtime: `{summary.get('runtime_s', 0):.2f}` s",
        f"- Cut rounds: `{summary.get('cut_rounds')}`",
        f"- Added cuts: `{summary.get('added_cuts')}`",
        f"- Input nodes/edges: `{summary.get('input_nodes')}` / `{summary.get('input_edges')}`",
        f"- Solution nodes/edges: `{summary.get('solution_nodes')}` / `{summary.get('solution_edges')}`",
        f"- Solution cycle rank: `{summary.get('solution_cycle_rank')}`",
        f"- Solution connected: `{summary.get('solution_connected')}`",
        f"- Solution bridge count: `{summary.get('solution_bridge_count')}`",
        f"- Solution is 2-edge-connected: `{summary.get('solution_is_2edge_connected')}`",
        f"- Source incident constraints: `{summary.get('source_connectivity_constraints')}`",
        f"- Source min incident edges: `{summary.get('source_min_incident_edges')}`",
        f"- Sources with insufficient candidates: `{summary.get('insufficient_source_degree_candidates')}`",
        f"- Original sources with candidate edges: `{summary.get('original_sources_with_candidate_edges')}`",
        f"- Selected original sources with incident edge: `{summary.get('selected_original_sources_with_incident_edge')}`",
        f"- Objective mode: `{summary.get('objective_mode')}`",
        f"- Coverage attribute: `{summary.get('coverage_attr')}`",
        f"- Chain extra cost budget: `{summary.get('chain_extra_cost_budget_m')}` m",
        f"- Selected chain extra cost: `{summary.get('selected_chain_extra_cost_m')}` m",
        f"- Total chain tree service length: `{summary.get('total_chain_tree_service_length_m')}` m",
        f"- Selected edge coverage: `{summary.get('selected_edge_coverage_weight', 0):.3f}`",
        f"- Total edge coverage available: `{summary.get('total_edge_coverage_weight', 0):.3f}`",
        f"- Selected edge coverage fraction: `{summary.get('selected_edge_coverage_fraction')}`",
        f"- Selected total coverage fraction: `{summary.get('selected_total_coverage_fraction')}`",
        f"- Objective length: `{summary.get('objective_length_m', 0):.3f}` m",
        f"- Best bound from phase `{summary.get('best_bound_phase')}`: `{best_bound_text}`",
        f"- MIP gap: `{summary.get('mip_gap')}`",
        f"- Redundancy constraint: `{summary.get('redundancy_constraint')}`",
        f"- Max redundancy constraint: `{summary.get('max_redundancy_constraint')}`",
        f"- Cut mode: `{summary.get('cut_mode')}`",
        f"- Threads: `{summary.get('threads')}`",
        f"- Requested MIPGap: `{summary.get('requested_mip_gap')}`",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Solve a minimum-length 2-edge-connected ILP on P2U candidates.")
    parser.add_argument("--input-gpkg", type=Path, default=INPUT_GPKG)
    parser.add_argument("--redundancy", type=int, default=None)
    parser.add_argument("--max-redundancy", type=int, default=None)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--mip-gap", type=float, default=0.02)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--max-cut-rounds", type=int, default=100)
    parser.add_argument("--cut-mode", choices=["iterative", "callback"], default="iterative")
    parser.add_argument(
        "--objective-mode",
        choices=["min_length", "max_chain_demand_then_min_length", "max_chain_demand_under_cost"],
        default="min_length",
    )
    parser.add_argument("--coverage-attr", default="edge_size_kva")
    parser.add_argument("--coverage-tolerance", type=float, default=1e-6)
    parser.add_argument("--source-min-incident-edges", type=int, default=1)
    parser.add_argument("--chain-extra-cost-budget-m", type=float, default=None)
    parser.add_argument("--output-gpkg", type=Path, default=OUTPUT_GPKG)
    parser.add_argument("--summary-json", type=Path, default=SUMMARY_JSON)
    parser.add_argument("--summary-md", type=Path, default=SUMMARY_MD)
    args = parser.parse_args()

    graph, edges, transformer_nodes, source_nodes, original_source_incident = load_ilp_graph(args.input_gpkg)
    solution, summary = solve_min_2edge(
        graph,
        original_source_incident=original_source_incident,
        redundancy=args.redundancy,
        max_redundancy=args.max_redundancy,
        time_limit=args.time_limit,
        mip_gap=args.mip_gap,
        threads=args.threads,
        max_cut_rounds=args.max_cut_rounds,
        cut_mode=args.cut_mode,
        objective_mode=args.objective_mode,
        coverage_attr=args.coverage_attr,
        coverage_tolerance=args.coverage_tolerance,
        source_min_incident_edges=args.source_min_incident_edges,
        chain_extra_cost_budget_m=args.chain_extra_cost_budget_m,
    )
    write_solution(
        solution,
        summary,
        edges,
        transformer_nodes,
        source_nodes,
        output_gpkg=args.output_gpkg,
        summary_json=args.summary_json,
        summary_md=args.summary_md,
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
