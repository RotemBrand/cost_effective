"""
Original-style real-data comparison with synthetic random switches.

This figure keeps the saved paper topologies from data/real_networks.nxjson:
original network versus optimized network. It then adds switch metadata as a
separate layer, using the same total switch counts observed in the SMART-DS SFO
feeders. The plotted indexes remain the original ones: Z_w, Z_R, and Z_F.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

if __package__ is None or __package__ == "":
    sys.path.append(str(Path(__file__).resolve().parents[1]))

import utilities.read_write as rw
from figures.real_data import DATA_PATH, get_sfo_bounded_network
from figures.real_data_plot import net_plot as original_net_plot
from indexes import Float, GraphRel, contract_switch_sections, edge_probs_by_length
from indexes.switch_placement import (
    add_synthetic_switches,
    count_switch_edges,
    count_tie_edges,
    ensure_edge_lengths,
)
from utilities.figures_utilities import BLUE, GREEN, GREY, RED, cm_to_inch


OUTPUT_DIR = Path("outputs") / "real_data_ties"
OUTPUT_PNG = OUTPUT_DIR / "sfo_random_switches_original_style.png"
OUTPUT_SVG = OUTPUT_DIR / "sfo_random_switches_original_style.svg"
OUTPUT_MD = OUTPUT_DIR / "sfo_random_switches_original_style_stats.md"

SFO_REFERENCE = {
    "SFO Davidson": ("P3U", "box"),
    "SFO Pacific": ("P2U", "box2"),
}

MEAN_EDGE_FAILURE_PROB = 5e-4
SIMULATION_SEED = 100
SWITCH_SEED = 2026
T_DAYS = 365 * 5
MEAN_CYCLE_DAYS = 0.01

TEXT_SIZE = 13
LABEL_SIZE = 12


def _edge_key(edge: tuple) -> tuple:
    return tuple(sorted(edge[:2]))


def _total_length(graph: nx.Graph) -> float:
    ensure_edge_lengths(graph)
    return sum(float(data["length_m"]) for _, _, data in graph.edges(data=True))


def _set_euclidean_lengths(graph: nx.Graph) -> None:
    pos = nx.get_node_attributes(graph, "pos")
    for u, v, data in graph.edges(data=True):
        if u in pos and v in pos:
            length = float(np.linalg.norm(np.asarray(pos[u], dtype=float) - np.asarray(pos[v], dtype=float)))
        else:
            length = float(data.get("length_m", data.get("length", data.get("weight", 1.0))))
        data["length_m"] = length
        data["length"] = length


def _reference_switch_count(name: str) -> int:
    area, box = SFO_REFERENCE[name]
    graph = get_sfo_bounded_network(area, box, include_ties=True)
    return count_switch_edges(graph)


def _power_rows() -> pd.DataFrame:
    networks_df = rw.read_nxjson(DATA_PATH)
    return networks_df[networks_df["type"] == "power"].copy()


def _prepare_graph(graph: nx.Graph, *, n_switches: int, seed: int) -> nx.Graph:
    graph = graph.copy()
    _set_euclidean_lengths(graph)
    if "sources" not in graph.graph:
        graph.graph["sources"] = [next(iter(graph.nodes))]
    return add_synthetic_switches(graph, n_switches=n_switches, seed=seed)


def _switch_aware_saidi(
    graph: nx.Graph,
    *,
    edge_failure_rate: float | None,
    rng: np.random.Generator,
) -> tuple[float, float]:
    sources = graph.graph.get("sources", [next(iter(graph.nodes))])
    section_graph = contract_switch_sections(
        graph,
        sources=sources,
        node_weight_attr="weight",
        edge_weight_attr="edge_weight",
        length_attr="length_m",
    )
    if edge_failure_rate is None:
        edge_probs, edge_failure_rate = edge_probs_by_length(
            section_graph,
            p=MEAN_EDGE_FAILURE_PROB,
            mode="mean",
            length_attr="length_m",
        )
    else:
        edge_probs, _ = edge_probs_by_length(
            section_graph,
            p=edge_failure_rate,
            mode="rate",
            length_attr="length_m",
        )
    edge_probs = {edge: Float(prob) for edge, prob in edge_probs.items()}
    graph_rel = GraphRel(section_graph, edges_prob=edge_probs, sources=sources)
    result = graph_rel.calc_rel_simulation(
        rel_type="saidi",
        T_days=T_DAYS,
        mean_cycle_days=MEAN_CYCLE_DAYS,
        rng=rng,
        show_progress=False,
    )
    return float(result.rel_result), float(edge_failure_rate)


def _stats_for_pair(row: pd.Series, rng: np.random.Generator) -> dict:
    name = row["name"]
    n_switches = _reference_switch_count(name)
    original = _prepare_graph(row["graph"], n_switches=n_switches, seed=SWITCH_SEED)
    optimized = _prepare_graph(row["optimal_network"], n_switches=n_switches, seed=SWITCH_SEED + 1)

    original_saidi, rate = _switch_aware_saidi(original, edge_failure_rate=None, rng=rng)
    optimized_saidi, _ = _switch_aware_saidi(optimized, edge_failure_rate=rate, rng=rng)

    original_weight = float(row["total_weight"])
    optimized_weight = float(row["optimal_network_weight"])
    original_r = original.number_of_edges() - original.number_of_nodes() + 1
    optimized_r = optimized.number_of_edges() - optimized.number_of_nodes() + 1

    return {
        "name": name,
        "original": original,
        "optimized": optimized,
        "reference_switches": n_switches,
        "original_switches": count_switch_edges(original),
        "optimized_switches": count_switch_edges(optimized),
        "original_ties": count_tie_edges(original),
        "optimized_ties": count_tie_edges(optimized),
        "N": original.number_of_nodes(),
        "M": original.number_of_edges(),
        "R": original_r,
        "R_optimized": optimized_r,
        "total_weight": original_weight,
        "optimal_network_weight": optimized_weight,
        "saidi": original_saidi,
        "optimal_network_saidi": optimized_saidi,
        "weight_ratio": 1 - optimized_weight / original_weight,
        "R_ratio": 1 - optimized_r / original_r if original_r else 0.0,
        "rel_ratio": 1 - optimized_saidi / original_saidi if original_saidi else 0.0,
        "edge_failure_rate": rate,
    }


def build_random_switch_comparison() -> pd.DataFrame:
    rng = np.random.default_rng(SIMULATION_SEED)
    stats = [_stats_for_pair(row, rng) for _, row in _power_rows().iterrows()]
    return pd.DataFrame(stats)


def smart_format(x: float) -> str:
    if abs(x) < 2:
        return f"{x:.2f}"
    if abs(x) < 10:
        return f"{x:.1f}"
    return f"{x:.0f}"


def _ratio_barplots(row: pd.Series, ax: plt.Axes) -> None:
    ratios = row[["weight_ratio", "R_ratio", "rel_ratio"]].values.astype(float)
    max_abs = max(np.max(np.abs(ratios)), 1e-12)
    norms = ratios / max_abs
    colors = [BLUE if ratio >= 0 else RED for ratio in ratios]
    labels = [r"$Z_w$", r"$Z_R$", r"$Z_F$"]
    y_pos = [2, 1, 0]

    ax.barh(y_pos, norms, color=colors, height=0.5)
    for val, norm, color, y in list(zip(ratios, norms, colors, y_pos))[::-1]:
        x = norm + 0.08 if norm >= 0 else norm - 0.08
        ha = "left" if norm >= 0 else "right"
        ax.text(x, y, smart_format(val), fontsize=TEXT_SIZE, color=color, va="center", ha=ha)

    ax.set_xlim((-1.25, 1.25))
    ax.set_ylim((-0.5, 2.5))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=TEXT_SIZE)
    ax.set_xticks([])
    ax.vlines(0, ymin=-0.25, ymax=2.25, linestyle="-", color=GREY)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _edge_midpoints(graph: nx.Graph, edgelist) -> tuple[list[float], list[float]]:
    pos = nx.get_node_attributes(graph, "pos")
    xs, ys = [], []
    for u, v in edgelist:
        if u not in pos or v not in pos:
            continue
        x = (pos[u][0] + pos[v][0]) / 2
        y = (pos[u][1] + pos[v][1]) / 2
        xs.append(x)
        ys.append(y)
    return xs, ys


def _net_plot(graph: nx.Graph, ax: plt.Axes) -> None:
    original_net_plot(graph, ax=ax, node_size=3, to_3857=False)

    closed_switch_edges = [
        (u, v)
        for u, v, data in graph.edges(data=True)
        if data.get("is_switch", False) and not data.get("is_tie", False)
    ]
    tie_edges = [(u, v) for u, v, data in graph.edges(data=True) if data.get("is_tie", False)]

    xs, ys = _edge_midpoints(graph, closed_switch_edges)
    if xs:
        ax.scatter(xs, ys, s=5, color=GREEN, alpha=0.55, linewidth=0, zorder=110)
    xs, ys = _edge_midpoints(graph, tie_edges)
    if xs:
        ax.scatter(xs, ys, s=18, color=BLUE, edgecolor="white", linewidth=0.25, zorder=120)


def _add_network_text(ax: plt.Axes, graph: nx.Graph, *, show_n: bool) -> None:
    n = graph.number_of_nodes()
    r = graph.number_of_edges() - graph.number_of_nodes() + 1
    if show_n:
        ax.text(0.05, 0.98, fr"$N = {{{n}}}$", ha="left", va="top", transform=ax.transAxes, fontsize=LABEL_SIZE)
    ax.text(0.95, 0.98, fr"$R = {{{r}}}$", ha="right", va="top", transform=ax.transAxes, fontsize=LABEL_SIZE)


def random_switches_original_style_plot(save: bool = True) -> plt.Figure:
    stats_df = build_random_switch_comparison()
    fig, axs = plt.subplots(
        len(stats_df),
        3,
        figsize=(3 * 8 / cm_to_inch, len(stats_df) * 6 / cm_to_inch),
        gridspec_kw={"wspace": 0.08},
        constrained_layout=True,
    )
    if len(stats_df) == 1:
        axs = axs.reshape(1, 3)

    for row_idx, (_, row) in enumerate(stats_df.iterrows()):
        _net_plot(row.original, axs[row_idx, 0])
        _net_plot(row.optimized, axs[row_idx, 1])
        _ratio_barplots(row, axs[row_idx, 2])

        axs[row_idx, 0].set_ylabel(row["name"], fontsize=TEXT_SIZE, color=RED)
        _add_network_text(axs[row_idx, 0], row.original, show_n=True)
        _add_network_text(axs[row_idx, 1], row.optimized, show_n=False)
        axs[row_idx, 0].text(
            0.05,
            0.08,
            f"switches={row.original_switches}, ties={row.original_ties}",
            transform=axs[row_idx, 0].transAxes,
            fontsize=LABEL_SIZE - 2,
        )
        axs[row_idx, 1].text(
            0.05,
            0.08,
            f"switches={row.optimized_switches}, ties={row.optimized_ties}",
            transform=axs[row_idx, 1].transAxes,
            fontsize=LABEL_SIZE - 2,
        )

    for ax, title in zip(axs[0, :2], ["Original", "Optimized"]):
        ax.text(0.5, 1.08, title, ha="center", va="bottom", transform=ax.transAxes, fontsize=TEXT_SIZE)

    handles = [
        plt.Line2D([0], [0], color=GREY, lw=1.2, label="Closed line"),
        plt.Line2D([0], [0], color=BLUE, lw=2.4, label="Tie/redundant edge"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=GREEN, markersize=5, label="Closed switch"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=RED, markeredgecolor="black", markersize=5, label="Source"),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.43, -0.03),
        ncol=4,
        frameon=False,
        fontsize=LABEL_SIZE - 1,
    )

    if save:
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
        fig.savefig(OUTPUT_SVG, bbox_inches="tight", transparent=True, format="svg")
        _write_markdown(stats_df)
    return fig


def _write_markdown(stats_df: pd.DataFrame) -> None:
    columns = [
        "name",
        "N",
        "M",
        "R",
        "R_optimized",
        "reference_switches",
        "original_switches",
        "optimized_switches",
        "original_ties",
        "optimized_ties",
        "total_weight",
        "optimal_network_weight",
        "saidi",
        "optimal_network_saidi",
        "weight_ratio",
        "R_ratio",
        "rel_ratio",
        "edge_failure_rate",
    ]
    lines = [
        "# SFO Random-Switch Original-Style Figure Stats",
        "",
        "Topology source: saved paper networks from `data/real_networks.nxjson`.",
        "",
        "Switch layer: non-tree edges are marked as tie switches; remaining switches are placed randomly on closed edges until the total switch count matches the corresponding SMART-DS feeder.",
        "",
        f"MCMC settings: `p={MEAN_EDGE_FAILURE_PROB:g}` mean section-edge failure probability for each original switched graph, `T_days={T_DAYS}`, `mean_cycle_days={MEAN_CYCLE_DAYS}`, simulation seed `{SIMULATION_SEED}`, switch seed `{SWITCH_SEED}`.",
        "",
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for _, row in stats_df.iterrows():
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, float):
                values.append(f"{value:.6g}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    lines.extend([
        "",
        "Index definitions:",
        "",
        "- `Z_w = weight_ratio = 1 - W_optimized / W_original`.",
        "- `Z_R = R_ratio = 1 - R_optimized / R_original`.",
        "- `Z_F = rel_ratio = 1 - F_optimized / F_original`.",
        "- `F` is computed after random switch placement and switch-section contraction.",
    ])
    OUTPUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    random_switches_original_style_plot(save=True)
