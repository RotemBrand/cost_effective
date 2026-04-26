import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon
import networkx as nx
from typing import List, Tuple, Dict
import optimal_network as ON
import optimal_network.other_methods.improve_tree as IT
import pandas as pd
from scipy.spatial.distance import pdist, squareform
import utilities.read_write as rw
from utilities.figures_utilities import (
    BLUE, RED, GREY, GREEN, YELLOW, cm_to_inch, TEXT_SIZE, LABEL_SIZE,
    LINE_WIDTH, configure_axes, set_custom_scientific_format, format_p_scientific,
    saidi_with_lengths, failing_rate_from_spanning_tree, save_axs_without_text
)

from indexes.graph_rel import GraphRel
from utilities.helper import draw_network
from indexes.utilities import create_sparse_graph
from indexes.utilities import get_skeleton_graph
from indexes.probs import Poly
from tqdm import tqdm

NODES_LINEWIDTH = 1.5
STRC_SIZE = 100
WIDTH = 1.5
NODES_EDGECOLOR='black'
WSPACE = 0.02


###### plot ######

def optimal_algorithm_plot(
        simulation_path: str = r"data/optimize_plot.nxjson",
        simulation_opt_path: str = r"data/optimize_plot_opt.nxjson",
        save: bool = False,
    ) -> None:
    """Generate the full optimization comparison figure."""
    import matplotlib as mpl

    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial"]
    # read data
    data = rw.read_nxjson(simulation_path).set_index(["category", "index"])

    # figures
    fig, axs = plt.subplots(
        5, 1,
        figsize=(3 * 8 /cm_to_inch, 5 * 8 / cm_to_inch),
        gridspec_kw={'hspace': 0.45}
    )
    rows = []

    ## illustration ##
    ax_main = axs[0]
    ax_row = make_ax_row(ax_main, 3)
    rows.append(ax_row)
    number_plots(ax_row, 0)
    xlim, ylim = draw_illustration(ax_row[0])
    draw_illustration_strc(ax_row[1])#,xlim=xlim, ylim=ylim)
    draw_chains_illustrations(ax_row[2])

    ## construction method ##
    ax_row = make_ax_row(axs[1], 3, margin=0.05, wspace=0.03)
    rows.append(ax_row)
    number_plots(ax_row, 3)
    with_colors = False
    draw_building_method(draw_clusters=True, draw_strc=False, draw_graph=False, with_colors=with_colors, ax=ax_row[0])
    draw_building_method(draw_clusters=True, draw_strc=True, draw_graph=False, with_colors=with_colors, ax=ax_row[1])
    draw_building_method(draw_clusters=True, draw_strc=False, draw_graph=True, with_colors=with_colors, ax=ax_row[2])
    titles = ["I. Clustering", r"II. Constructing $G_{struct}$", r"III. Constructing $u_q$"]
    for ax, title in zip(ax_row, titles):
        ax.axis(False)
        # ax.text(0.5, 1.07, title, ha='center', va='bottom', color=BLUE, fontsize=TEXT_SIZE)


    ## trees ##
    ax_row = make_ax_row(axs[2], 3)
    rows.append(ax_row)
    number_plots(ax_row, 6)
    plot_tree_row(ax_row, data)

    ## non optimal ##
    ax_row = make_ax_row(axs[3], 4, wspace=[WSPACE, WSPACE, WSPACE*6])
    rows.append(ax_row)
    number_plots(ax_row, 9)
    plot_non_optimal_row(ax_row, data)

    ## optimal ##
    data_opt = rw.read_nxjson(simulation_opt_path)
    ax_main = axs[4]
    ax_row = make_ax_row(ax_main, 3, wspace=[WSPACE, WSPACE*6])
    rows.append(ax_row)
    number_plots(ax_row, 12)
    plot_optimal_row(ax_row, data, data_opt)

    # text
    add_text(rows)

    # save
    if save:
        fig.savefig(
            r'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\optimize_plot\optimize_plot.svg',
            bbox_inches='tight',
            dpi=300,
            transparent=True
        )
        axs_row = [ax for row in rows for ax in row]
        save_axs_without_text(axs_row, r'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\optimize_plot\optimize_plot\axes')

# iluustrations
def draw_illustration(ax: plt.Axes) -> tuple[tuple, tuple]:
    """Draw the structure graph illustration."""
    graph = create_illustrative_graph()
    comps = list(nx.k_edge_components(graph, 2))
    node_to_color = {
        node: YELLOW if node == 0
        else RED if len(comp) == 1
        else BLUE if graph.degree(node) >= 3
        else GREEN
        for comp in comps
        for node in comp
    }
    node_to_size = {
        node: 15 * 1.2 if len(comp) == 1
        else 60 * 1.2 if graph.degree(node) >= 3
        else 35 * 1.2
        for comp in comps
        for node in comp
    }
    draw_network(
        graph,
        node_color=[YELLOW if node == 0 else 'white' for node in graph.nodes],
        node_size=[node_to_size[node] for node in graph.nodes],
        edge_color='grey',
        # node_size=[200 if node == 0 else 100 for node in graph.nodes],
        with_labels=True,
        labels={0: r'$s$'},
        font_color='black',
        width=WIDTH,
        alpha=1,
        linewidths=NODES_LINEWIDTH,
        edgecolors=[node_to_color[node] for node in graph.nodes],
        ax=ax
    )
    add_legend_to_illustration(ax)
    return ax.get_xlim(), ax.get_ylim()

def add_legend_to_illustration(ax: plt.Axes) -> None:
    """Add a legend for trees, chains, and forks to the illustration."""
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=RED,   lw=3, label='Trees'),
        Line2D([0], [0], color=GREEN, lw=3, label='Chains'),
        Line2D([0], [0], color=BLUE,  lw=3, label='Forks'),
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper left',
        bbox_to_anchor=(0.0, 1.2), # (x,y) anchor point in axes coords
        bbox_transform=ax.transAxes,  # explicit: use axes [0..1] coords
        fontsize=TEXT_SIZE - 2.5,
        frameon=False,
        handlelength=1.0,     # length of the sample line
        labelcolor='linecolor',
    )

def _key(e: tuple) -> tuple:
    if len(e) == 2:
        return tuple(sorted(e))
    return (*sorted(e[:2]), e[2])

def draw_illustration_strc(ax: plt.Axes, xlim: tuple = None) -> None:
    """Draw the skeleton structure for the illustration graph."""
    # get graph structure
    graph = create_illustrative_graph()
    comps = list(nx.k_edge_components(graph, 2))
    nodes = set([node for comp in comps if len(comp) > 1 for node in comp])
    subgraph = nx.subgraph(graph, nodes).copy()
    strc = get_skeleton_graph(subgraph)
 
    # width of each chain is propto its width
    edge_to_chain_length = {
        _key(e): data["length"]
        for chain, data in strc.edges.items()
        for e in data["subgraph"].edges
    }
    widths = np.array([edge_to_chain_length[_key(e)] / 1.5 for e in subgraph.edges])

    # draw
    nx.draw_networkx_edges(
        subgraph,
        pos=nx.get_node_attributes(subgraph, "pos"),
        width=widths / widths.max() * WIDTH * 2,
        edge_color=GREEN,
        ax=ax

    )
    nx.draw_networkx_nodes(
        subgraph,
        pos=nx.get_node_attributes(subgraph, "pos"),
        node_size=[0 if d == 2 else 120 for node, d in subgraph.degree],
        node_color=[YELLOW if node == 0 else 'white' for node in subgraph.nodes],
        edgecolors=[YELLOW if node == 0 else BLUE for node in subgraph.nodes],
        linewidths=NODES_LINEWIDTH,
        ax=ax
    )
    nx.draw_networkx_labels(
        subgraph,
        pos=nx.get_node_attributes(subgraph, "pos"),
        labels={0: r'$s$'},
        font_color='black',
        ax=ax
    )

    # configure
    ax.axis(False)
    if xlim:
        ax.set_xlim(xlim)
        ax.set_ylim(xlim)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin - (ymax - ymin) * 0.1, ymax)


def draw_chains_illustrations(ax: plt.Axes) -> None:
    """Draw the individual chains for the illustration."""
    lengths = [6, 3, 6, 1]
    max_length = max(lengths)
    xmax = max_length + 2
    # center = xmax / 2
    for i, length in enumerate(lengths):
        # xleft = center - length  / 2 - 1
        # xright = center + length / 2 - 1
        # x = np.linspace(xleft, xright, length +2)
        row = len(lengths) - i
        x = np.arange(1.25, length + 3, 1)
        y = [row] * len(x)
        print(len(x))
        ax.plot(
            x,
            y,
            color=GREY,
            linewidth=WIDTH,
            zorder=1,
        )
        sizes = [0] + [100] * (len(x) - 2) + [0]
        ax.scatter(
            x,
            y,
            s=sizes,
            facecolors='white',
            edgecolors=GREEN,
            linewidths=NODES_LINEWIDTH,
            zorder=2,
        )
        ax.text(0, row, fr'$u_{{{i}}}$', fontsize=TEXT_SIZE, color=GREEN, ha='left')
        # configure ax
        ax.axis(False)
        ax.set_ylim((0, len(lengths) + 1))

# Explain the method
def draw_building_method(draw_clusters: bool, draw_strc: bool, draw_graph: bool, with_colors: bool, ax: plt.Axes) -> None:
    """Draw the steps for constructing the optimal network (clustering, strc, full graph)."""
    r = 4
    n = 100
    seed=42

    # create data
    rng = np.random.default_rng(seed)
    points = rng.random((n, 2))

    n = 2 * (r - 1)
    m = 3 * (r - 1)

    # clusters the points into chains
    chains, centers = ON.balanced_kmeans_gurobi(
        points, m, max_iter=7, chain_len_sigma=0, random_state=seed
    )

    data = pd.DataFrame(np.c_[points, chains], columns=["x", "y", "c"])
    data['node'] = list(range(len(points)))
    data['node'] = data['node'].apply(int).astype('int')
    data['c'] = data['c'].astype(int)

    # pick an ordered palette and bind each cluster id -> color
    clusters = list(set(data['c']))
    if with_colors:
        color_map = dict(zip(clusters, sns.color_palette(n_colors=len(clusters))))
    else:
        color_map = {c: GREEN for c in clusters}
    data['color'] = data['c'].map(color_map)

    # scatter with the exact same palette mapping (legend hidden as before)
    if draw_clusters:
        if not draw_graph:
            sns.scatterplot(
                x=data['x'], y=data['y'], c='white', edgecolors='GREEN',
                s=30, ax=ax, zorder=1, linewidths=NODES_LINEWIDTH
            )

        # add convex hull for each cluster
        for c in clusters:
            pts = data.loc[data['c'] == c, ['x', 'y']].to_numpy()
            if len(pts) < 3:
                continue  # need at least 3 points for a hull
            hull = ConvexHull(pts)
            poly_xy = pts[hull.vertices]
            ax.add_patch(Polygon(
                poly_xy, closed=True,
                facecolor=color_map[c],
                alpha=0.20,
                edgecolor=color_map[c],
                linewidth=1.5,
                zorder=0
            ))


    ######## trips ########
    chosen_trips = ON.get_optimal_strc_trips(centers, points, n_nodes=n, strc_n_init_iters=2, exact_vertices=False, debug=False, source_node=None)

    trips_G = nx.Graph()
    for s, t, c in chosen_trips:
        cs = str(c)
        trips_G.add_edge(s, cs)
        trips_G.add_edge(cs, t)
        trips_G.nodes[s]["pos"] = points[s]
        trips_G.nodes[t]["pos"] = points[t]
        trips_G.nodes[cs]["pos"] = centers[c]
    if draw_strc:
        draw_network(
            trips_G,
            node_color=['white' if d != 2 else GREY for node, d in trips_G.degree],
            edgecolors=[GREY if d == 2 else BLUE for node, d in trips_G.degree],
            edge_color='grey',
            linewidths=NODES_LINEWIDTH,
            chain_size=20,
            strc_size=STRC_SIZE,
            with_labels=False,
            width=WIDTH + 1,
            ax=ax,
        )
    
    ######## graph ########
    if draw_graph:
        subgraph = ON.add_chains_to_strc(points, chains, chosen_trips, max_init_iter=2, debug=False)
        
        draw_network(
            subgraph,
            node_color='white',
            # node_size=0,
            strc_size=STRC_SIZE,
            chain_size=30,
            with_labels=False,
            edge_color='grey',
            linewidths=NODES_LINEWIDTH,
            edgecolors=[GREEN if d == 2 else BLUE for node, d in subgraph.degree],
            width=WIDTH,
            ax=ax
        )


# examples
def draw_tree_graph(tree: nx.Graph, ax: plt.Axes) -> None:
    """Draw a tree graph."""
    draw_network(
        graph=tree,
        edge_color=[RED if d['source'] == "tree" else BLUE for e, d in tree.edges.items()],
        node_size=0,
        with_labels=False,
        width=WIDTH,
        ax=ax,
    )

def draw_2con_optimal(graph: nx.Graph, ax: plt.Axes, width: float = WIDTH, chain_color: str = GREEN, arc_edges: bool = False) -> None:
    """Draw a 2-connected optimal graph or structure."""
    draw_network(
        graph=graph,
        node_color='white',
        edgecolors=[BLUE if d >=  3 else chain_color for node, d in graph.degree],
        node_size=[50 if d >=  3 else 0 for node, d in graph.degree],
        edge_color=chain_color,
        with_labels=False,
        linewidths=NODES_LINEWIDTH,
        width=width,
        arc_edges=arc_edges,
        ax=ax
    )


def draw_fr_line(data_opt: pd.DataFrame, ax: plt.Axes) -> None:
    """Draw the F vs R line for the optimal graphs."""
    r = data_opt["r"].values.astype('float')
    f = data_opt["saidi"].values.astype('float')
    sns.lineplot(x=[0, max(r)], y=0, ax=ax, color='black', linewidth=LINE_WIDTH, zorder=0)
    sns.scatterplot(
        x=r, y=f, ax=ax, s=60, zorder=100,
        # color=BLUE,
        facecolor='white',
        edgecolor=BLUE,
        linewidth=1.5,
    )
    sns.lineplot(
        x=r, y= f[0] * r ** (-2), ax=ax, linewidth=LINE_WIDTH
    )
    ax.text(0, 0, r"$F_{3Reg}$", ha='left', va='bottom', fontsize=TEXT_SIZE-2, zorder=200, color="black")
    ax.set_xlabel(r"$R$", fontsize=LABEL_SIZE, color=BLUE)
    ax.set_ylabel(r"$F_{opt}$", fontsize=LABEL_SIZE, color=BLUE)
    ax.tick_params(axis='both', colors=BLUE)
    configure_axes(ax)

def add_arrow(ax: plt.Axes, x: float, y: float, l: float) -> FancyArrowPatch:
    """Add a horizontal arrow to the plot."""
    arrow = FancyArrowPatch(
        (x, y), (x + l, y),
        transform=ax.transAxes,
        arrowstyle='-|>',           # line with a filled head
        color='black',
        linewidth=LINE_WIDTH,
        mutation_scale=0.5,
        # shrinkA=0, shrinkB=0,
        clip_on=False,
        # zorder=zorder,
        # **arrowprops
    )
    ax.add_patch(arrow)
    return arrow

# plot trees
def plot_tree_row(ax_row: list[plt.Axes], data: pd.DataFrame) -> None:
    """Plot the row of tree graphs."""

    ## trees ##
    original_graph = data.loc[("graph", 0)]["graph"]
    original_graph.remove_edges_from(original_graph.edges)
    draw_network(
        original_graph,
        node_color=GREY,
        node_size=5,
        with_labels=False,
        edgecolors=GREY,
        ax=ax_row[0]
    )
    tree_graph, tree_saidi, tree_p = data.loc[("trees", 0)][["graph", "saidi", "p"]]
    draw_tree_graph(tree_graph, ax=ax_row[1])


    ax = ax_row[2]
    values = [tree_saidi, tree_p]
    bars = ax.bar([r'$F$', r'$p$'], values, color=RED, width=0.5)
    
    # add value text on top of bars
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2, bar.get_height(),
            format_p_scientific(val),
            ha="center", va="bottom",
            color=RED, fontsize=TEXT_SIZE,
        )
    # format ax
    ax.set_yticks([])
    ax.set_ylim((0, max(values) * 1.3))
    ax.tick_params(axis='x', colors=RED)
    
    for spine in ax.spines.values():
        spine.set_visible(False)

def plot_non_optimal_row(ax_row: list[plt.Axes], data: pd.DataFrame) -> None:
    """Plot the row of the non optimal graphs."""
    # plot graphs
    non_optimal = data.loc[("non_optimals", slice(None))]
    non_optimal["sigma"] = non_optimal["graph"].apply(_get_effective_sigma)
    for i, (graph, saidi, p, sigma) in enumerate(non_optimal[["graph", "saidi", "p", "sigma"]].values):
        ax = ax_row[i]
        draw_2con_optimal(graph, ax=ax, chain_color=GREEN)
        ax.text(
            0.5, -0.02, rf"$\sigma_{{\tilde{{\lambda}}}}$={sigma:.1f}",
            ha='center', va='top', fontsize=TEXT_SIZE, color=GREEN
        )

    # horizontal bar plot
    ax = ax_row[3]


    sigmas = non_optimal["sigma"].values
    values = non_optimal["saidi"].values
    y = np.arange(len(sigmas))

    bars = ax.barh(y, values, color=GREEN, height=0.8)

    # add value text to the right of bars
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_width(),                          # x-position
            bar.get_y() + bar.get_height() / 2,       # y-position
            format_p_scientific(val),
            ha="left", va="center",
            color=GREEN, fontsize=TEXT_SIZE,
        )

    # format ax
    ax.set_xticks([])
    ax.set_yticks(y, [rf"${sigma:.0f}$" for sigma in sigmas])
    ax.set_xlim(0, max(values) * 1.3)
    ax.set_xlabel(r"$F$")
    ax.set_ylabel(r"$\sigma_{\tilde{\lambda}}$")
    for spine in ["right", "top"]:
        ax.spines[spine].set_visible(False)
    

    configure_axes(ax)

def plot_optimal_row(ax_row: list[plt.Axes], data: pd.DataFrame, data_opt: pd.DataFrame) -> None:
    """Plot the row of the optimal graphs."""
    optimal_row = data.loc[("optimal", 0)]

    # plot graph
    draw_2con_optimal(optimal_row["graph"], ax_row[0], chain_color=BLUE)
    sigma = _get_effective_sigma(optimal_row["graph"])
    ax_row[0].text(
        0.5, -0.02, rf"$\sigma_{{\tilde{{\lambda}}}}$={sigma:.1f}",
        ha='center', va='top', fontsize=TEXT_SIZE, color=BLUE
    )

    # plot strc
    strc = get_skeleton_graph(optimal_row["graph"], sources=[node for node, d in optimal_row["graph"].degree if d == 3])
    draw_2con_optimal(strc, ax_row[1], width=WIDTH+1, chain_color=BLUE, arc_edges=True)

    # optimal f vs f
    draw_fr_line(data_opt, ax_row[2])

def _get_effective_sigma(graph: nx.Graph) -> float:
    """Calculate the effective standard deviation of the chain lengths."""
    # get lengths
    for e in graph.edges:
        graph.edges[e]["length"] = np.linalg.norm(
            np.array(graph.nodes[e[0]]["pos"]) -
            np.array(graph.nodes[e[1]]["pos"])
        )
    # get strc
    strc = get_skeleton_graph(graph)
    lambdas = []
    for e, e_data in strc.edges.items():
        subgraph = e_data["subgraph"]
        total_length = sum([graph.edges[l]["length"] for l in subgraph.edges])
        total_weight = len(subgraph.nodes) - 2
        lambdas.append(total_length * np.sqrt(total_weight))
    return np.std(lambdas)


# plot utils

def make_ax_row(
    ax_main: plt.Axes,
    n_axs: int,
    margin: float = 0.0,
    wspace: float | list[float] = WSPACE,
) -> list[plt.Axes]:
    """
    Create a row of n_axs inset axes inside ax_main.

    Parameters
    ----------
    ax_main : matplotlib axis
        The parent axis.
    n_axs : int
        Number of child axes.
    margin : float, default=0.0
        Left/right margin as fraction of parent width.
    wspace : float or list of float, default=0.02
        Horizontal space(s) between axes as fraction of parent width.
        - If float: same spacing used between all axes.
        - If list: must have length n_axs-1, giving spacing between consecutive axes.

    Returns
    -------
    axs : list of matplotlib axes
        List of inset axes aligned in a row.
    """
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    # normalize wspace to a list
    if isinstance(wspace, (int, float)):
        wspaces = [float(wspace)] * (n_axs - 1)
    else:
        if len(wspace) != n_axs - 1:
            raise ValueError("If wspace is a list, it must have length n_axs - 1")
        wspaces = list(map(float, wspace))

    total_wspace = sum(wspaces)
    total_space = 1 - 2 * margin - total_wspace
    if total_space <= 0:
        raise ValueError("Margins + spaces exceed figure width")

    width_each = total_space / n_axs

    axs = []
    x0 = margin
    for i in range(n_axs):
        ax = inset_axes(
            ax_main,
            width="100%",
            height="100%",
            loc="lower left",
            bbox_to_anchor=(x0, 0, width_each, 1),
            bbox_transform=ax_main.transAxes,
            borderpad=0,
        )
        axs.append(ax)

        if i < n_axs - 1:
            x0 += width_each + wspaces[i]

    ax_main.axis(False)
    return axs

LETTERS = ['a', 'b', 'c' ,'d', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o']
def number_plots(ax_row: list[plt.Axes], start: int = 0) -> None:
    """Add letter labels (a, b, c...) to a row of subplots."""
    for i, ax in enumerate(ax_row):
        ax.text(
            -0.03, 1.05, f"{LETTERS[start + i]}",
            fontsize=TEXT_SIZE, color="black",
            ha='right', va='bottom',
            transform=ax.transAxes,
        )

def plot_line_with_title(ax: plt.Axes, title: str, color: str) -> None:
    """Plot a vertical line with a title text on the given axis."""
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x = xlim[0] - (xlim[1] - xlim[0]) *  0.07
    y = (ylim[1] + ylim[0]) / 2
    x_text = x - (xlim[1] - xlim[0]) * 0.01
    ax.axvline(
        [x], ymin=ylim[0]*0.95 + ylim[1]*0.05, ymax=ylim[0]*0.05 + ylim[1]*0.95,
        color=color, clip_on=False, linewidth=LINE_WIDTH
    )
    ax.text(
        x_text, y, title, rotation=90,
        ha='right', va='center', color=color, fontsize=TEXT_SIZE
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

###### simulations ######

def simulate_netowrks(file_path: str) -> tuple[pd.DataFrame, np.array, float]:
    """Simulate reliability for random trees, non-optimal, and optimal networks."""
    n = 1000
    r = 10
    c = int((n - 2 * (r - 1)) / (3 * (r - 1)))
    print(f"c = {c}")
    seed=42
    rng = np.random.default_rng(seed)
    points, graph = random_points(n=n, seed=seed)

    # simulate tree
    p_rate = failing_rate_from_spanning_tree(points, p_mean=5e-4)
    
    tree = nx.minimum_spanning_tree(graph, weight="weight")
    subgraph, _ = IT.improve_tree(tree, R=r+1, root=0, weight_quantile=1)
    
    # Tag edges for plotting
    for e in subgraph.edges:
        subgraph.edges[e]["source"] = "tree" if e in tree.edges else "chord"
        
    trees = [[subgraph, "trees", 0, get_saidi(subgraph, p_rate=p_rate, rng=rng)]]
    print("Simulate tree")

    # simulate non-optimal 2-connected
    non_optimals = []
    sigmas = [c // 4, c, 4 * c]
    for i, sigma in enumerate(sigmas):
        subgraph = ON.optimal_network_from_points(
            points=points,
            r=r,
            kmeans_max_iter=7,
            chain_len_sigma=sigma,
            seed=7,
            strc_n_init_iters=2,
        )
        non_optimals.append([subgraph, "non_optimals", sigma, get_saidi(subgraph, p_rate=p_rate, rng=rng)])
        print(f"Simulate non-optimal ({i+1}/3)")

    # simulate optimal
    optimal_subgraph = ON.optimal_network_from_points(
        points=points,
        r=r,
        kmeans_max_iter=7,
        chain_len_sigma=0,
        seed=7,
        strc_n_init_iters=2,
    )
    optimal = [[optimal_subgraph, "optimal", 0, get_saidi(optimal_subgraph, p_rate=p_rate, rng=rng)]]
    print(f"Simulate optimal")
    
    # original_graph
    original_graph = [[graph, "graph", 0, None]]
    data = pd.DataFrame(
        original_graph + trees + non_optimals + optimal,
        columns=["graph", "category", "index", "saidi"]
    )
    data['p'] = [5e-4 for _ in range(len(data))]
    data['p_rate'] = [p_rate for _ in range(len(data))]
    # save
    rw.write_nxjson(data, file_path)
    return data, points, p_rate

def simulate_opt_networks(p_rate: float, points: np.array, file_name: str=r"data/graph_five_opt.nxjson") -> pd.DataFrame:
    """Simualte optimal network for each r in range(r_min, r_max, step) and return a df of r and saidi"""
    n = 1000
    r_list = range(1, 11, 1)
    rng = np.random.default_rng(2000)
    
    # calculate saidi for each optimal r
    data = [] # ["graph", "source", "saidi", "r", "p", "p_rate"]
    for r in tqdm(r_list):
        # build graph
        optimal_subgraph = ON.optimal_network_from_points(
            points=points,
            r=r,
            kmeans_max_iter=7,
            chain_len_sigma=0,
            seed=7,
            strc_n_init_iters=2,
        )
        # saidi
        saidi = get_saidi(optimal_subgraph, p_rate=p_rate, rng=rng)
        data.append([optimal_subgraph.copy(), 0, saidi, r, 5e-4, p_rate])
    
    # save
    data = pd.DataFrame(data, columns=["graph", "source", "saidi", "r", "p", "p_rate"])
    rw.write_nxjson(data, file_name)
    return data



def get_saidi(G: nx.Graph, p_rate: float, rng: np.random.Generator) -> list:
    # source = next(iter([node for node, d in G.degree if d == 3]))
    source = 0
    saidi = saidi_with_lengths(G, sources=[source], p=p_rate, mode="rate", rng=rng, mean_cycle_days=0.5)
    return saidi

def random_points(n: int, seed: int = None) -> tuple[np.array, nx.Graph]:
    """Generate random points and build a base proximity graph."""
    # create data
    rng = np.random.default_rng(seed)
    points = rng.random((n, 2))
    dist_mat = squareform(pdist(points))
    graph = nx.from_numpy_array(dist_mat)
    edges_to_remove = [(u, v) for u, v, weight in graph.edges(data='weight') if weight > 0.2]
    graph.remove_edges_from(edges_to_remove)
    # pos
    pos = dict(enumerate(points))
    nx.set_node_attributes(graph, pos, "pos")
    return points, graph

def create_illustrative_graph() -> nx.Graph:
    """create the network to illustrate the structure graph"""
    # build strc
    strc = nx.from_edgelist([
        (0, 1), (0, 1), (1, 2), (2, 3), (2, 3), (3, 4), (4, 5), (5, 0), (0, 3)
    ], create_using=nx.MultiGraph)
    # add chains
    G = nx.Graph(create_sparse_graph(
        strc,
        {
            (0, 1, 0): 0, (0, 1, 1): 3, (1, 2, 0): 1,
            (2, 3, 0): 1, (2, 3, 1): 6, (3, 4, 0): 2,
            (4, 5, 0): 1, (5, 0, 0): 1,
        }
    ))
    # add trees
    G.add_edges_from([
        (2, 200), (2, 201),
        (4, 400),
        (5, 500), (500, 501), (500, 502), (502, 503),
        (500, 504), (504, 505), (504, 506), (506, 507), (507, 508), (508, 509)

    ])
    # pos
    pos = nx.kamada_kawai_layout(G, center=(0, 0))
    pos = {
        node: (-y, x)
        for node, (x, y) in pos.items()
    }
    def _add(node: int, x: float, y: float):
        p = pos[node]
        pos[node] = (p[0] + x, p[1] + y)

    _add(20, 0.1, 0)
    _add(502, 0.05, .1)
    _add(503, 0, .25)
    _add(501, .15, 0)
    _add(505, 0, .1)
    _add(200, -0.1, 0.2)
    _add(201, -0.1, -0.2)
    nx.set_node_attributes(G, pos, "pos")
    return G


##### text ###
def add_text(axs: list[list[plt.Axes]]) -> None:
    """Add text annotations to the axes."""
    # b
    ax = axs[0][1]
    ax.text(
        0, 1.02, r"$G_{struct}$", transform=ax.transAxes,
        ha='left', va="top", fontsize=TEXT_SIZE, color=BLUE
    )
    ax.text(
        0.5, 0, r"$F_{struct}$", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE, color=BLUE
    )
    ax.text(
        1, 0, r"$\ll$", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE, color=BLUE
    )

    # c
    ax = axs[0][2]
    ax.text(
        0.5, 1.02, "Intermediate chains", transform=ax.transAxes,
        ha='center', va="top", fontsize=TEXT_SIZE, color=BLUE
    )
    ax.text(
        0.5, 0, r"$F_{inter}$", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE, color=BLUE
    )

    # d e f
    ax_row = axs[1]
    texts = ["I. Clustering", r"II. Constructing $G_{strc}$", r"III. Constructing $u_{q}$"]
    for ax, text in zip(ax_row, texts):
        ax.text(
            0.5, 1.05, text, transform=ax.transAxes,
            ha='center', va="bottom", fontsize=TEXT_SIZE - 2.5, color=BLUE
        )
    
    # g
    ax = axs[2][0]
    ax.text(
        0.5, 1.05, r"$N=1000, R=10$", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE - 2.5, color=GREY
    )

    # h
    ax = axs[2][1]
    ax.text(
        0.5, 1.05, "Naive", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE, color=RED
    )

    # k
    ax = axs[3][1]
    ax.text(
        0.5, 1.2, "2-connected", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE, color=GREEN
    )

    # m
    ax = axs[4][0]
    ax.text(
        1, 1.2, "Optimal", transform=ax.transAxes,
        ha='center', va="bottom", fontsize=TEXT_SIZE, color=BLUE
    )

    # n
    ax = axs[4][1]
    ax.text(
        0.1, 0.95, r"$G_{struct}$", transform=ax.transAxes,
        ha='left', va="bottom", fontsize=TEXT_SIZE, color=BLUE
    )

        
