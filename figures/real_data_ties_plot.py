"""
Side figure for inspecting SMART-DS SFO feeders with normally open tie switches.

This script intentionally does not read or overwrite data/real_networks.nxjson.
It rebuilds the two SFO networks from the raw SMART-DS files and saves a
separate tie-aware diagnostic figure under outputs/real_data_ties/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(str(Path(__file__).resolve().parents[1]))

from figures.real_data import (
    get_ext_chains,
    get_optimal_saidi,
    get_prec_2_conn,
    get_sfo_bounded_network,
)
from indexes import Float, GraphRel, defult_sources, edge_probs_by_length
from optimal_network import optimal_network_from_points
from utilities.figures_utilities import BLUE, GREEN, GREY, RED, cm_to_inch


OUTPUT_DIR = Path("outputs") / "real_data_ties"
OUTPUT_PNG = OUTPUT_DIR / "sfo_tie_switches.png"
OUTPUT_SVG = OUTPUT_DIR / "sfo_tie_switches.svg"
OUTPUT_MD = OUTPUT_DIR / "sfo_tie_switches_stats.md"
OPT_OUTPUT_PNG = OUTPUT_DIR / "sfo_tie_switches_optimized.png"
OPT_OUTPUT_SVG = OUTPUT_DIR / "sfo_tie_switches_optimized.svg"
OPT_OUTPUT_MD = OUTPUT_DIR / "sfo_tie_switches_optimized_stats.md"

SFO_NETWORKS = (
    ("SFO Davidson", "P3U", "box"),
    ("SFO Pacific", "P2U", "box2"),
)

TEXT_SIZE = 10
STATS_SIZE = 7
TABLE_SIZE = 6.5
MEAN_EDGE_FAILURE_PROB = 5e-4
SIMULATION_SEED = 100
OPTIMIZATION_SEED = 10
T_DAYS = 365 * 5
MEAN_CYCLE_DAYS = 0.01


def _edge_key(edge: tuple) -> tuple:
    return tuple(sorted(edge[:2]))


def _edge_lists(graph: nx.Graph) -> tuple[list[tuple], list[tuple], list[tuple]]:
    tie_edges = [
        (u, v)
        for u, v, data in graph.edges(data=True)
        if data.get("is_tie", False)
    ]
    switch_edges = [
        (u, v)
        for u, v, data in graph.edges(data=True)
        if data.get("is_switch", False) and not data.get("is_tie", False)
    ]
    regular_edges = [
        (u, v)
        for u, v, data in graph.edges(data=True)
        if not data.get("is_tie", False)
    ]
    return regular_edges, switch_edges, tie_edges


def _draw_edges(
    graph: nx.Graph,
    edgelist: list[tuple],
    ax: plt.Axes,
    *,
    color,
    width: float,
    alpha: float = 1.0,
    linestyle: str = "solid",
    zorder: int = 1,
):
    from matplotlib.collections import LineCollection

    pos = nx.get_node_attributes(graph, "pos")
    segments = [
        [(pos[u][0], pos[u][1]), (pos[v][0], pos[v][1])]
        for u, v in edgelist
        if u in pos and v in pos
    ]
    if not segments:
        return None

    collection = LineCollection(
        segments,
        colors=[color] * len(segments),
        linewidths=width,
        linestyles=linestyle,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_collection(collection)
    return collection


def _set_limits(graphs: list[nx.Graph], axs):
    xs, ys = [], []
    for graph in graphs:
        for x, y in nx.get_node_attributes(graph, "pos").values():
            xs.append(x)
            ys.append(y)
    if not xs or not ys:
        return
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    x_pad = (x_max - x_min) * 0.04
    y_pad = (y_max - y_min) * 0.04
    for ax in axs:
        ax.set_xlim(x_min - x_pad, x_max + x_pad)
        ax.set_ylim(y_min - y_pad, y_max + y_pad)


def _draw_sources(graph: nx.Graph, ax: plt.Axes):
    pos = nx.get_node_attributes(graph, "pos")
    sources = [node for node in graph.graph.get("sources", []) if node in pos]
    if not sources:
        return
    ax.scatter(
        [pos[node][0] for node in sources],
        [pos[node][1] for node in sources],
        s=22,
        color=RED,
        edgecolor="black",
        linewidth=0.4,
        zorder=5,
    )


def _format_axes(ax: plt.Axes):
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def _count_ties(graph: nx.Graph) -> int:
    return sum(1 for _, _, data in graph.edges(data=True) if data.get("is_tie", False))


def _count_switches(graph: nx.Graph) -> int:
    return sum(1 for _, _, data in graph.edges(data=True) if data.get("is_switch", False))


def _total_length_m(graph: nx.Graph) -> float:
    return sum(float(data.get("length_m", data.get("length", 0.0))) for _, _, data in graph.edges(data=True))


def _edge_length_from_pos(graph: nx.Graph, u, v) -> float:
    pos = nx.get_node_attributes(graph, "pos")
    return float(np.linalg.norm(np.asarray(pos[u], dtype=float) - np.asarray(pos[v], dtype=float)))


def _ensure_edge_lengths(graph: nx.Graph) -> None:
    for u, v, data in graph.edges(data=True):
        if "length_m" not in data:
            data["length_m"] = _edge_length_from_pos(graph, u, v)
        data.setdefault("length", data["length_m"])


def _mark_redundant_edges_as_ties(graph: nx.Graph) -> None:
    """Mark the non-tree edges of an abstract optimized graph as tie switches."""
    _ensure_edge_lengths(graph)
    tree_edges = {_edge_key(edge) for edge in nx.minimum_spanning_tree(graph, weight="length_m").edges}
    for u, v, data in graph.edges(data=True):
        is_tie = _edge_key((u, v)) not in tree_edges
        data["is_tie"] = is_tie
        data["is_switch"] = is_tie
        data["normally_closed"] = not is_tie


def _source_index(graph: nx.Graph) -> int | None:
    sources = list(graph.graph.get("sources", []))
    if not sources:
        return None
    nodes = list(graph.nodes)
    return nodes.index(sources[0])


def _optimized_same_r_graph(graph: nx.Graph, *, seed: int) -> nx.Graph:
    """Build an optimized graph on the same node positions and with the same R."""
    redundancy = graph.number_of_edges() - graph.number_of_nodes() + 1
    if redundancy < 1:
        raise ValueError("The optimized comparison requires a graph with R >= 1")

    nodes = list(graph.nodes)
    points = np.asarray([graph.nodes[node]["pos"] for node in nodes], dtype=float)
    source_idx = _source_index(graph)
    optimized = optimal_network_from_points(
        points,
        r=redundancy,
        source_node=source_idx,
        strc_n_init_iters=3,
        chain_n_init_iters=3,
        seed=seed,
        debug=False,
    )
    optimized.graph["sources"] = [source_idx] if source_idx is not None else []
    _mark_redundant_edges_as_ties(optimized)
    return optimized


def _network_stats(graph: nx.Graph) -> str:
    tie_count = sum(1 for _, _, data in graph.edges(data=True) if data.get("is_tie", False))
    switch_count = sum(1 for _, _, data in graph.edges(data=True) if data.get("is_switch", False))
    return (
        f"N={graph.number_of_nodes()}, M={graph.number_of_edges()}, "
        f"R={graph.number_of_edges() - graph.number_of_nodes() + 1}, "
        f"ties={tie_count}, switches={switch_count}"
    )


def _simulate_saidi(graph: nx.Graph, edge_failure_rate: float, rng: np.random.Generator) -> float:
    edge_probs, _ = edge_probs_by_length(
        graph,
        p=edge_failure_rate,
        mode="rate",
        length_attr="length_m",
    )
    edge_probs = {edge: Float(prob) for edge, prob in edge_probs.items()}
    sources = graph.graph.get("sources", [defult_sources(graph)])
    graph_rel = GraphRel(graph, edges_prob=edge_probs, sources=sources)
    result = graph_rel.calc_rel_simulation(
        rel_type="saidi",
        T_days=T_DAYS,
        mean_cycle_days=MEAN_CYCLE_DAYS,
        rng=rng,
        show_progress=False,
    )
    return float(result.rel_result)


def _stats_for_graph(
    network_name: str,
    case: str,
    graph: nx.Graph,
    *,
    edge_failure_rate: float,
    rng: np.random.Generator,
) -> dict:
    n_nodes = graph.number_of_nodes()
    n_edges = graph.number_of_edges()
    redundancy = n_edges - n_nodes + 1
    chains = get_ext_chains(nx.Graph(graph))
    total_length_m = _total_length_m(graph)
    n_sources = len(graph.graph.get("sources", []))
    return {
        "network": network_name,
        "case": case,
        "N": n_nodes,
        "M": n_edges,
        "R": redundancy,
        "R/N": redundancy / n_nodes if n_nodes else 0.0,
        "sources": n_sources,
        "ties": _count_ties(graph),
        "switches": _count_switches(graph),
        "total_length_km": total_length_m / 1000,
        "mean_edge_length_m": total_length_m / n_edges if n_edges else 0.0,
        "optimal_c": (n_nodes - 2 * (redundancy - 1)) / (3 * (redundancy - 1)) if redundancy > 1 else n_nodes,
        "optimal_saidi": get_optimal_saidi(n_nodes, redundancy, MEAN_EDGE_FAILURE_PROB, n_sources or 1),
        "prec2_conn": get_prec_2_conn(nx.Graph(graph)),
        "ext_chains_len": chains,
        "max_chain": max(chains) if chains else 0,
        "top_chains": chains[:5],
        "saidi": _simulate_saidi(graph, edge_failure_rate, rng),
    }


def _build_network_rows() -> tuple[list[tuple[str, nx.Graph, nx.Graph]], list[dict]]:
    rng = np.random.default_rng(SIMULATION_SEED)
    graph_rows = []
    stats_rows = []

    for name, area, box in SFO_NETWORKS:
        radial = get_sfo_bounded_network(area, box, include_ties=False)
        with_ties = get_sfo_bounded_network(area, box, include_ties=True)
        graph_rows.append((name, radial, with_ties))

        _, edge_failure_rate = edge_probs_by_length(
            radial,
            p=MEAN_EDGE_FAILURE_PROB,
            mode="mean",
            length_attr="length_m",
        )
        stats_rows.append(
            _stats_for_graph(
                name,
                "Radial feeder",
                radial,
                edge_failure_rate=edge_failure_rate,
                rng=rng,
            )
        )
        stats_rows.append(
            _stats_for_graph(
                name,
                "With normally open ties",
                with_ties,
                edge_failure_rate=edge_failure_rate,
                rng=rng,
            )
        )

    return graph_rows, stats_rows


def _format_value(value: float, *, digits: int = 3) -> str:
    if isinstance(value, int):
        return str(value)
    if abs(value) >= 100:
        return f"{value:.1f}"
    if abs(value) >= 10:
        return f"{value:.2f}"
    return f"{value:.{digits}g}"


def _stats_by_network(stats_rows: list[dict]) -> dict[str, dict[str, dict]]:
    by_network: dict[str, dict[str, dict]] = {}
    for row in stats_rows:
        by_network.setdefault(row["network"], {})[row["case"]] = row
    return by_network


def _comparison_gains(original_stats: dict, optimized_stats: dict) -> dict:
    return {
        "Z_w": 1 - optimized_stats["total_length_km"] / original_stats["total_length_km"],
        "Z_R": 1 - optimized_stats["R"] / original_stats["R"] if original_stats["R"] else 0.0,
        "Z_F": 1 - optimized_stats["saidi"] / original_stats["saidi"] if original_stats["saidi"] else 0.0,
    }


def _add_stats_table(ax: plt.Axes, radial_stats: dict, tie_stats: dict) -> None:
    ax.axis("off")
    metrics = [
        ("N", "N"),
        ("M", "M"),
        ("R", "R"),
        ("R/N", "R/N"),
        ("ties", "ties"),
        ("switches", "switches"),
        ("length km", "total_length_km"),
        ("2-conn", "prec2_conn"),
        ("c*", "optimal_c"),
        ("F sim", "saidi"),
        ("F*/ideal", "optimal_saidi"),
        ("max chain", "max_chain"),
    ]
    cell_text = [
        [
            label,
            _format_value(radial_stats[key]),
            _format_value(tie_stats[key]),
        ]
        for label, key in metrics
    ]
    table = ax.table(
        cellText=cell_text,
        colLabels=["metric", "radial", "ties"],
        cellLoc="center",
        colLoc="center",
        bbox=[0.0, 0.03, 1.0, 0.9],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(TABLE_SIZE)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor((0.85, 0.85, 0.85))
        if row == 0:
            cell.set_facecolor((0.94, 0.94, 0.94))
            cell.set_text_props(weight="bold")


def _add_optimized_stats_table(ax: plt.Axes, original_stats: dict, optimized_stats: dict) -> None:
    ax.axis("off")
    gains = _comparison_gains(original_stats, optimized_stats)
    metrics = [
        ("N", "N", ""),
        ("M", "M", ""),
        ("R", "R", "Z_R"),
        ("R/N", "R/N", ""),
        ("ties", "ties", ""),
        ("length km", "total_length_km", "Z_w"),
        ("2-conn", "prec2_conn", ""),
        ("c*", "optimal_c", ""),
        ("F sim", "saidi", "Z_F"),
        ("F*/ideal", "optimal_saidi", ""),
        ("max chain", "max_chain", ""),
    ]
    cell_text = []
    for label, key, gain_key in metrics:
        gain = _format_value(gains[gain_key]) if gain_key else ""
        cell_text.append([
            label,
            _format_value(original_stats[key]),
            _format_value(optimized_stats[key]),
            gain,
        ])
    table = ax.table(
        cellText=cell_text,
        colLabels=["metric", "orig.", "opt.", "gain"],
        cellLoc="center",
        colLoc="center",
        bbox=[0.0, 0.03, 1.0, 0.9],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(TABLE_SIZE)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor((0.85, 0.85, 0.85))
        if row == 0:
            cell.set_facecolor((0.94, 0.94, 0.94))
            cell.set_text_props(weight="bold")


def _markdown_float(value: float, digits: int = 6) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.{digits}g}"


def _write_stats_markdown(stats_rows: list[dict]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    columns = [
        "network",
        "case",
        "N",
        "M",
        "R",
        "R/N",
        "sources",
        "ties",
        "switches",
        "total_length_km",
        "mean_edge_length_m",
        "optimal_c",
        "optimal_saidi",
        "prec2_conn",
        "saidi",
        "max_chain",
        "top_chains",
    ]
    lines = [
        "# SFO Tie-Switch Figure Stats",
        "",
        "Generated by `figures/real_data_ties_plot.py` from SMART-DS SFO feeder data.",
        "",
        f"MCMC settings: mean edge failure probability target `{MEAN_EDGE_FAILURE_PROB:g}` on the radial feeder, "
        f"`T_days={T_DAYS}`, `mean_cycle_days={MEAN_CYCLE_DAYS}`, seed `{SIMULATION_SEED}`. "
        "The tie-enabled case uses the same length-based failure rate as its radial counterpart.",
        "",
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in stats_rows:
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, float):
                values.append(_markdown_float(value))
            elif isinstance(value, list):
                values.append(", ".join(map(str, value)))
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    lines.extend([
        "",
        "Metric notes:",
        "",
        "- `R = M - N + 1`.",
        "- `total_length_km` is the sum of SMART-DS `length_m` edge attributes.",
        "- `prec2_conn` is the fraction of nodes in nontrivial 2-edge-connected components.",
        "- `optimal_c` and `optimal_saidi` use the same helper definitions as `figures.real_data.networks_dict_to_df`.",
        "- `saidi` is the MCMC disconnected load fraction using default unit node weights.",
    ])
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_optimized_markdown(stats_rows: list[dict]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    columns = [
        "network",
        "case",
        "N",
        "M",
        "R",
        "R/N",
        "sources",
        "ties",
        "switches",
        "total_length_km",
        "mean_edge_length_m",
        "optimal_c",
        "optimal_saidi",
        "prec2_conn",
        "saidi",
        "max_chain",
        "top_chains",
    ]
    lines = [
        "# SFO Original-vs-Optimized Tie-Switch Stats",
        "",
        "Generated by `figures/real_data_ties_plot.py` from SMART-DS SFO feeder data.",
        "",
        "The optimized graph is built on the same node positions as the original tie-enabled feeder "
        "and uses the same redundancy index `R = M - N + 1`. For the optimized abstract graph, "
        "non-tree edges are marked as tie switches for plotting and counting.",
        "",
        f"MCMC settings: mean edge failure probability target `{MEAN_EDGE_FAILURE_PROB:g}` on the original feeder, "
        f"`T_days={T_DAYS}`, `mean_cycle_days={MEAN_CYCLE_DAYS}`, simulation seed `{SIMULATION_SEED}`, "
        f"optimization seed `{OPTIMIZATION_SEED}`. The optimized case uses the same length-based failure rate as its original feeder.",
        "",
        "## Network Stats",
        "",
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in stats_rows:
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, float):
                values.append(_markdown_float(value))
            elif isinstance(value, list):
                values.append(", ".join(map(str, value)))
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")

    lines.extend(["", "## Original-to-Optimized Gains", "", "| network | Z_w | Z_R | Z_F |", "| --- | --- | --- | --- |"])
    by_network = _stats_by_network(stats_rows)
    for network, rows in by_network.items():
        gains = _comparison_gains(rows["Original with ties"], rows["Optimized same R"])
        lines.append(
            "| "
            + " | ".join([
                network,
                _markdown_float(gains["Z_w"]),
                _markdown_float(gains["Z_R"]),
                _markdown_float(gains["Z_F"]),
            ])
            + " |"
        )

    lines.extend([
        "",
        "Metric notes:",
        "",
        "- `Z_w = 1 - length_optimized / length_original`.",
        "- `Z_R = 1 - R_optimized / R_original`; this should be zero here because `R` is fixed.",
        "- `Z_F = 1 - F_optimized / F_original`, where `F` is the MCMC disconnected load fraction.",
        "- `total_length_km` is the sum of edge `length_m`; optimized edge lengths are Euclidean lengths from node positions.",
        "- `prec2_conn` is the fraction of nodes in nontrivial 2-edge-connected components.",
        "- `optimal_c` and `optimal_saidi` use the same helper definitions as `figures.real_data.networks_dict_to_df`.",
        "- `saidi` is the MCMC disconnected load fraction using default unit node weights.",
    ])
    OPT_OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _plot_network(graph: nx.Graph, ax: plt.Axes, *, show_switches: bool):
    regular_edges, switch_edges, tie_edges = _edge_lists(graph)
    switch_edge_keys = {_edge_key(switch_edge) for switch_edge in switch_edges}
    non_switch_regular = [
        edge for edge in regular_edges
        if _edge_key(edge) not in switch_edge_keys
    ]

    _draw_edges(graph, non_switch_regular, ax, color=GREY, width=0.7, alpha=0.55, zorder=1)
    if show_switches:
        _draw_edges(graph, switch_edges, ax, color=GREEN, width=0.9, alpha=0.65, zorder=2)
    else:
        _draw_edges(graph, switch_edges, ax, color=GREY, width=0.7, alpha=0.55, zorder=1)
    _draw_edges(graph, tie_edges, ax, color=BLUE, width=2.4, alpha=1.0, linestyle="dashed", zorder=4)
    _draw_sources(graph, ax)
    _format_axes(ax)


def _build_optimized_rows() -> tuple[list[tuple[str, nx.Graph, nx.Graph]], list[dict]]:
    rng = np.random.default_rng(SIMULATION_SEED)
    graph_rows = []
    stats_rows = []

    for network_idx, (name, area, box) in enumerate(SFO_NETWORKS):
        original = get_sfo_bounded_network(area, box, include_ties=True)
        optimized = _optimized_same_r_graph(original, seed=OPTIMIZATION_SEED + network_idx)
        graph_rows.append((name, original, optimized))

        _, edge_failure_rate = edge_probs_by_length(
            original,
            p=MEAN_EDGE_FAILURE_PROB,
            mode="mean",
            length_attr="length_m",
        )
        stats_rows.append(
            _stats_for_graph(
                name,
                "Original with ties",
                original,
                edge_failure_rate=edge_failure_rate,
                rng=rng,
            )
        )
        stats_rows.append(
            _stats_for_graph(
                name,
                "Optimized same R",
                optimized,
                edge_failure_rate=edge_failure_rate,
                rng=rng,
            )
        )

    return graph_rows, stats_rows


def sfo_tie_switch_plot(save: bool = True) -> plt.Figure:
    """Create the side-by-side SFO tie-switch diagnostic figure."""
    rows, stats_rows = _build_network_rows()
    stats_by_network = _stats_by_network(stats_rows)

    fig, axs = plt.subplots(
        len(rows),
        3,
        figsize=(24 / cm_to_inch, 15 / cm_to_inch),
        gridspec_kw={"width_ratios": [1.0, 1.0, 0.62]},
        constrained_layout=True,
    )
    if len(rows) == 1:
        axs = axs.reshape(1, 3)

    for row_index, (name, radial, with_ties) in enumerate(rows):
        _set_limits([radial, with_ties], axs[row_index, :2])

        _plot_network(radial, axs[row_index, 0], show_switches=False)
        _plot_network(with_ties, axs[row_index, 1], show_switches=True)
        _add_stats_table(
            axs[row_index, 2],
            stats_by_network[name]["Radial feeder"],
            stats_by_network[name]["With normally open ties"],
        )

        axs[row_index, 0].set_ylabel(name, fontsize=TEXT_SIZE, labelpad=8)
        for ax, graph in zip(axs[row_index, :], (radial, with_ties)):
            ax.text(
                0.02,
                0.98,
                _network_stats(graph),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=STATS_SIZE,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
            )

    axs[0, 0].set_title("Radial feeder", fontsize=TEXT_SIZE + 1)
    axs[0, 1].set_title("With normally open ties", fontsize=TEXT_SIZE + 1)
    axs[0, 2].text(0.5, 1.01, "Stats", transform=axs[0, 2].transAxes, ha="center", fontsize=TEXT_SIZE + 1)

    handles = [
        plt.Line2D([0], [0], color=GREY, lw=1.2, label="Closed line"),
        plt.Line2D([0], [0], color=GREEN, lw=1.2, label="Switchable closed line"),
        plt.Line2D([0], [0], color=BLUE, lw=2.4, linestyle="dashed", label="Normally open tie"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markeredgecolor="black", markersize=5, label="Source"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.42, -0.01),
        ncol=4,
        frameon=False,
        fontsize=TEXT_SIZE - 1,
    )

    if save:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
        fig.savefig(OUTPUT_SVG, bbox_inches="tight", transparent=True, format="svg")
        _write_stats_markdown(stats_rows)
    return fig


def sfo_tie_optimized_plot(save: bool = True) -> plt.Figure:
    """Create an original-with-ties vs optimized-same-R comparison figure."""
    rows, stats_rows = _build_optimized_rows()
    stats_by_network = _stats_by_network(stats_rows)

    fig, axs = plt.subplots(
        len(rows),
        3,
        figsize=(24 / cm_to_inch, 15 / cm_to_inch),
        gridspec_kw={"width_ratios": [1.0, 1.0, 0.62]},
        constrained_layout=True,
    )
    if len(rows) == 1:
        axs = axs.reshape(1, 3)

    for row_index, (name, original, optimized) in enumerate(rows):
        _set_limits([original, optimized], axs[row_index, :2])

        _plot_network(original, axs[row_index, 0], show_switches=True)
        _plot_network(optimized, axs[row_index, 1], show_switches=True)
        _add_optimized_stats_table(
            axs[row_index, 2],
            stats_by_network[name]["Original with ties"],
            stats_by_network[name]["Optimized same R"],
        )

        axs[row_index, 0].set_ylabel(name, fontsize=TEXT_SIZE, labelpad=8)
        for ax, graph in zip(axs[row_index, :2], (original, optimized)):
            ax.text(
                0.02,
                0.98,
                _network_stats(graph),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=STATS_SIZE,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 1.5},
            )

    axs[0, 0].set_title("Original with ties", fontsize=TEXT_SIZE + 1)
    axs[0, 1].set_title("Optimized same R", fontsize=TEXT_SIZE + 1)
    axs[0, 2].text(0.5, 1.01, "Stats", transform=axs[0, 2].transAxes, ha="center", fontsize=TEXT_SIZE + 1)

    handles = [
        plt.Line2D([0], [0], color=GREY, lw=1.2, label="Closed line"),
        plt.Line2D([0], [0], color=GREEN, lw=1.2, label="Switchable closed line"),
        plt.Line2D([0], [0], color=BLUE, lw=2.4, linestyle="dashed", label="Tie/redundant edge"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markeredgecolor="black", markersize=5, label="Source"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.42, -0.01),
        ncol=4,
        frameon=False,
        fontsize=TEXT_SIZE - 1,
    )

    if save:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(OPT_OUTPUT_PNG, dpi=300, bbox_inches="tight")
        fig.savefig(OPT_OUTPUT_SVG, bbox_inches="tight", transparent=True, format="svg")
        _write_optimized_markdown(stats_rows)
    return fig


if __name__ == "__main__":
    sfo_tie_optimized_plot(save=True)
