"""
This file is the main real-networks plot
based on the data generated in real_data.py and saves into real_data.DATA_PATH
"""
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from pathlib import Path
import contextily as ctx
import seaborn as sns
import networkx as nx
import pandas as pd
import numpy as np
import utilities.read_write as rw
from .real_data import DATA_PATH
from utilities.helper import draw_network
from utilities.figures_utilities import cm_to_inch, BLUE, GREY, RED, GREEN, _pos_4326_to_3857


TEXT_SIZE = 16
c_by_type = {
    'communication': GREEN,
    'water': BLUE,
    'power': RED
}

# load type icons (optional; falls back gracefully if images missing)
ICON_DIR = Path(__file__).parent
_type_icons = {}
ICON_ZOOM = 0.06  # adjust if you want larger/smaller icons
TEXT_OFFSET = 0.04  # horizontal offset for text when icon present
for k in c_by_type.keys():
    p = ICON_DIR / f"{k}.png"
    if p.exists():
        try:
            _type_icons[k] = mpimg.imread(str(p))
        except Exception:
            # ignore read errors and continue without icon
            pass



def real_data_plot(show_basemap: bool=True, show_text: bool=True,  save=False):
    """
    The main function for the real data plot
    """
    
    # read data
    networks_df = rw.read_nxjson(DATA_PATH)
    from IPython.display import display
    print(networks_df.columns)
    display(
        networks_df[[
            "name", "total_weight", "optimal_network_weight", "N", "M", "R"
            ,'R_ratio', 'weight_ratio', 'rel_ratio', 'optimal_network_saidi', 'saidi'
            ]]
        )

    # configure fig
    fig, axs = plt.subplots(6, 3, figsize=(3*8/cm_to_inch, 6*6/cm_to_inch))
    set_grid_spacing(fig, axs, odd_hspace=0.025, even_hspace=0.055, even_wspace=0.06, odd_wspace=0.04)

    # plot networks
    for i, (_, row) in enumerate(networks_df.iterrows()):
        color = c_by_type[row.type]
        net_plot(row.graph, ax=axs[i, 0], node_size=6 if row.type == "communication" else 3, to_3857=row.type == "communication")
        net_plot(row.optimal_network, ax=axs[i, 1], node_size=6 if row.type == "communication" else 3, to_3857=row.type == "communication")
        if show_basemap:
            if row.type == "communication":
                ctx.add_basemap(axs[i, 0], crs=3857, source=ctx.providers.CartoDB.Positron)
                ctx.add_basemap(axs[i, 1], crs=3857, source=ctx.providers.CartoDB.Positron)
            elif row.type == "power":
                ctx.add_basemap(axs[i, 0], crs=3857, source=ctx.providers.OpenStreetMap.Mapnik)
                ctx.add_basemap(axs[i, 1], crs=3857, source=ctx.providers.OpenStreetMap.Mapnik)
            elif row.type == "water":
                for ax in axs[i, :2]:
                    ax.set_facecolor("#dbd8d8")

        # N, R labels
        ax = axs[i, 0]
        ax.text(0.05, 0.98, fr"$N = {{{row.N}}}$", ha='left', va='top', transform=ax.transAxes, fontsize=TEXT_SIZE)
        ax.text(0.95, 0.98, fr"$R = {{{row.R}}}$", ha='right', va='top', transform=ax.transAxes, fontsize=TEXT_SIZE)
        ax = axs[i, 1]
        R = row.R_optimal_network
        ax.text(0.95, 0.98, fr"$R = {{{R}}}$", ha='right', va='top', transform=ax.transAxes, fontsize=TEXT_SIZE)

        # ylabel
        axs[i, 0].axis(True)
        axs[i, 1].axis(True)
        axs[i, 0].set_ylabel(row["name"] , fontsize=TEXT_SIZE, color=color)

        # bar plots
        ratio_barplots(row, axs[i, 2])

    # titles
    for ax, title in zip(axs[0, :], ["Original", "Optimized"]):
        ax.text(
            0.5, 1.2, title, ha='center', va='bottom', transform=ax.transAxes, fontsize=TEXT_SIZE
        )
        
        # type titles
        for t, ax in zip(networks_df["type"].iloc[::2], axs[::2, 0]):

            ax.text(
                0, 1.03, t.capitalize(), ha='left', va='bottom', transform=ax.transAxes, fontsize=TEXT_SIZE, color=c_by_type[t]
            )

    # remove text
    for ax in axs[:, :2].flatten():
        for text in ax.texts:
            if text._text.startswith("(C)"):
                text._text = ""
    # numbers
    number_plots(axs=axs, start=0)
    # save
    if save:
        fig.savefig(
            r'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\real_data\real_data.svg',
            bbox_inches='tight',
            transparent=True,
            format='svg',
        )



def smart_format(x):
    """
    Format number with adaptive decimal digits:
    - < 1 → 2 decimals (e.g. 0.97 → '0.97')
    - < 10 → 1 decimal (e.g. 1.97 → '2.0')
    - >= 10 → 0 decimals (e.g. 20.34 → '20')
    """
    if x < 2:
        return f"{x:.2f}"
    elif x < 10:
        return f"{x:.1f}"
    else:
        return f"{x:.0f}"


def ratio_barplots(row: pd.Series, ax: plt.Axes):
    """
    Plot that stats barplots on the side of the main figure
    """
    # data
    ratios = row[["weight_ratio", "R_ratio", "rel_ratio"]].values.astype(float)
    print(ratios)
    norms = ratios.copy()
    norms /= norms.max()
    # norms[2:] /= norms[2:].max()
    colors = [BLUE if r >= 0 else RED for r in ratios]
    labels = [r"$Z_w$", r"$Z_R$", r"$Z_F$"]

    # custom y positions with gap
    y_pos = [2, 1, 0]

    ax.barh(y_pos, norms, color=colors, height=0.5)

    for val, norm, color, y in list(zip(ratios, norms, colors, y_pos))[::-1]:
        x = norm + 0.1 if norm >= 0 else 0.05
        ha = 'left'
        ax.text(
            x, y, smart_format(val),
            fontsize=TEXT_SIZE, color=color,
            va='center', zorder=100, ha=ha,
        )
    # x axis

    # configure axes
    ax.set_xlim((-0.15, 1.5))
    ax.set_ylim((-0.5, 2.5))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=TEXT_SIZE)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.vlines(0, ymin=-0.25, ymax=2.25, linestyle="-", color=GREY)



def split_axis_horizontally(ax, gap_ratio=0.04):
    """
    Split an existing Axes (inside a larger figure) into two horizontal sub-axes.
    Each new axis covers half the width and full height of the original one.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The existing Axes to split.
    gap_ratio : float
        Fractional gap between left and right axes (relative to original width).

    Returns
    -------
    ax_left, ax_right : matplotlib.axes.Axes
        The two new sub-axes occupying the original space.
    """
    fig = ax.figure
    pos = ax.get_position()  # original [x0, y0, width, height] in figure coords

    # Compute geometry
    gap = pos.width * gap_ratio
    half_width = (pos.width - gap) / 2

    left_rect = [pos.x0, pos.y0, half_width, pos.height]
    right_rect = [pos.x0 + half_width + gap, pos.y0, half_width, pos.height]

    # Hide the original Axes
    ax.axis(False)

    # Add two new Axes inside the same region
    ax_left = fig.add_axes(left_rect)
    ax_right = fig.add_axes(right_rect)

    return ax_left, ax_right


def set_grid_spacing(
    fig,
    axs,
    *,
    odd_hspace=0.06,
    even_hspace=0.03,
    odd_wspace=0.04,
    even_wspace=0.02,
):
    """
    Adjust per-gap spacing in a grid of axes by shifting axes positions.

    Spacing is in *figure coordinates* (0..1), i.e. the actual gap between axes
    bounding boxes.

    Vertical gaps:
      gap between row i and i+1 uses odd_hspace if i is even (0->1,2->3,...) else even_hspace.

    Horizontal gaps:
      gap between col j and j+1 uses odd_wspace if j is even (0->1,2->3,...) else even_wspace.

    Notes:
    - Works best when all axes have equal sizes initially (created via plt.subplots).
    - This shifts axes; it does not resize them unless shifts accumulate near edges.
    """
    axs = np.asarray(axs)
    nrows, ncols = axs.shape

    # Reset global spacings first (optional but helps make behavior deterministic)
    fig.subplots_adjust(hspace=0, wspace=0)

    # ---------- VERTICAL (row gaps) ----------
    # We shift rows downward to increase gaps (or upward if decreasing).
    for i in range(nrows - 1):
        # gap between row i (below) and row i+1 (above in matplotlib coords)?:
        # In figure coords, higher y is up. For typical subplots:
        # row i is above row i+1. So use row i+1 top and row i bottom carefully.
        # Let's compute gap between row i (upper) and row i+1 (lower):
        bottom_upper = axs[i, 0].get_position().y0
        top_lower = axs[i + 1, 0].get_position().y1
        current_gap = bottom_upper - top_lower

        target_gap = odd_hspace if (i % 2 == 0) else even_hspace
        delta = target_gap - current_gap

        if abs(delta) < 1e-12:
            continue

        # To increase gap, move all rows below downward by delta
        # (moving down decreases y). If delta is negative, we move up.
        for r in range(i + 1, nrows):
            for c in range(ncols):
                ax = axs[r, c]
                pos = ax.get_position()
                ax.set_position([pos.x0, pos.y0 - delta, pos.width, pos.height])

    # ---------- HORIZONTAL (column gaps) ----------
    # We shift columns to the right to increase gaps.
    for j in range(ncols - 1):
        # gap between col j (left) and col j+1 (right):
        right_left = axs[0, j].get_position().x1
        left_right = axs[0, j + 1].get_position().x0
        current_gap = left_right - right_left

        target_gap = odd_wspace if (j % 2 == 0) else even_wspace
        delta = target_gap - current_gap

        if abs(delta) < 1e-12:
            continue

        # To increase gap, move all columns to the right by delta
        # (moving right increases x). If delta is negative, we move left.
        for c in range(j + 1, ncols):
            for r in range(nrows):
                ax = axs[r, c]
                pos = ax.get_position()
                ax.set_position([pos.x0 + delta, pos.y0, pos.width, pos.height])

LETTERS = ['a', 'b', 'c' ,'d', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r']
def number_plots(axs, start: int=0):
    for i, ax in enumerate(axs.flatten()):
        ax.text(
            -0.04, 1.0, f"{LETTERS[start + i]}",
            fontsize=TEXT_SIZE, color="black",
            ha='right', va='bottom',
            transform=ax.transAxes,
        )



# maps

def _key(e: tuple) -> tuple:
    return tuple(sorted(e))


def draw_edges(graph: nx.Graph, pos=None, edgelist=None, width=1.0, edge_color='k', ax=None, alpha=1.0, zorder=100, linestyle='solid'):
    """Draw edges using a LineCollection so we can control z-order and appearance.

    Parameters
    ----------
    graph : networkx.Graph
        Graph containing the edges (used only to infer nodes if edgelist is None).
    pos : dict or None
        Mapping node -> (x, y). If None, attempts to read node 'pos' attributes.
    edgelist : iterable of edge tuples or None
        List of edges to draw. If None, draws all edges in the graph.
    width : float or sequence
        Line width (can be scalar or per-edge list).
    edge_color : color or sequence
        Matplotlib color or list of colors per edge.
    ax : matplotlib.axes.Axes or None
        Axes to draw on. If None, uses current axes.
    alpha : float
        Alpha transparency for edges.
    zorder : int
        Drawing z-order (higher renders on top).
    linestyle : str
        Line style passed to LineCollection.
    """
    from matplotlib.collections import LineCollection
    if ax is None:
        ax = plt.gca()

    if pos is None:
        pos = nx.get_node_attributes(graph, 'pos')
    if edgelist is None:
        edgelist = list(graph.edges)

    # Build segments
    segments = []
    for u, v in edgelist:
        if u not in pos or v not in pos:
            continue
        p1 = pos[u]
        p2 = pos[v]
        segments.append([(p1[0], p1[1]), (p2[0], p2[1])])

    if len(segments) == 0:
        return None

    # Normalize colors and widths
    colors = edge_color
    lw = width

    lc = LineCollection(segments, colors=colors, linewidths=lw, linestyles=linestyle, alpha=alpha, zorder=zorder)
    ax.add_collection(lc)
    return lc

def net_plot(graph: nx.Graph, ax, node_size: float=3.0, to_3857: bool=False):
    sources = graph.graph.get("sources", set())
    T_edges = set(map(_key, nx.minimum_spanning_tree(graph).edges))
    pos = nx.get_node_attributes(graph, "pos")
    if to_3857:
        pos = _pos_4326_to_3857(pos)
    draw_network(
        graph,
        pos=pos,
        node_color=[RED if node in sources else 'black' for node in graph.nodes],
        node_size=[20 if node in sources else node_size for node in graph.nodes],
        edgecolors=None,
        with_labels=False,
        ax=ax
    )
    re_edges = [e for e in graph.edges if _key(e) not in T_edges]
    draw_edges(
        graph,
        pos=pos,
        edgelist=re_edges,
        width=3,
        edge_color=[BLUE] * len(re_edges),
        ax=ax,
        zorder=100,
    )

    ylim = np.array(ax.get_ylim())
    ylim[1] += (ylim[1] - ylim[0]) * 0.1
    ax.set_ylim(ylim)

        
