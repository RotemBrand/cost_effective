import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import networkx as nx
import utilities.read_write as rw
from utilities import draw_network
from utilities.figures_utilities import (
    RED, GREEN, GREY, BLUE, YELLOW,
    TEXT_SIZE, LABEL_SIZE, TICK_SIZE, ANNOTATION_SIZE,
    cm_to_inch, format_p_scientific, set_custom_scientific_format,
    configure_axes, net_plot, draw_edges
)

LINEWIDTH = 2.5
POINTSIZE=50

NEW_YORK_FILE_NAME = r'data/ny/manhattan.nxjson'
NEW_YORK_ASYMPTOTIC_FILE_NAME = r'data/ny/nyc_asymptotic.nxjson'
NEW_YORK_MCMC_FILE_NAME = r'data/ny/mcmc_simulation.nxjson'

def nature_ny_graph(
    new_york_file_name: str = NEW_YORK_FILE_NAME,
    new_york_asymptotic_file_name: str = NEW_YORK_ASYMPTOTIC_FILE_NAME,
    MCMC_file_name: str = NEW_YORK_MCMC_FILE_NAME,
    save: bool = False, 
) -> tuple[plt.Figure, dict]:
    """
    Generate the main New York network figure for Nature submission.
    
    Creates a 3-row layout figure showing tree and ring topologies, 
    reliability plots (MCMC and polynomial), and scaling analysis.
    """
    # figure
    fig, axs = make_figure_layout(figsize=(15, 17.5), group_gap=0.5, big_ratio=1.15, inter_group_gap=0.4, top_hspace=0.35)
    ny_data = rw.read_nxjson(new_york_file_name)


    # --- first row ---
    # tree network
    ax = axs['a']
    tree_row = ny_data.query('r == 0').iloc[0]
    from indexes import GraphRel
    rel = GraphRel.reliability_polynomial(tree_row.graph, max_fail=1, sources=tree_row.sources)
    net_plot(tree_row.graph, ax=ax, sources=tree_row.sources, basemap=True)

    # tree MCMCM
    ax = axs['b']
    reliability_MCMC_plot(data_file_name=MCMC_file_name, ax=ax, r=0)

    # tree p vs f
    ax = axs['c']
    p_vs_f_lineplot(
        data=ny_data, r_list=[0], ax=ax, color=BLUE, r_text_ha='left', xmax=1e-4
    )

    # ring network
    ax = axs['d']
    ring_row = ny_data.query('r == 1').iloc[0]
    net_plot(ring_row.graph, ax=ax, sources=ring_row.sources, basemap=True)

    # trringee MCMCM
    ax = axs['e']
    reliability_MCMC_plot(data_file_name=MCMC_file_name, ax=ax, r=1)

    # ring p vs f
    ax = axs['f']
    p_vs_f_lineplot(
        data=ny_data, r_list=[1], ax=ax, color=BLUE, r_text_ha='left', xmax=1e-4
    )

    # --- second row ---
    # mesh p vs f
    ax = axs['g']
    p_vs_f_lineplot(
        data=ny_data, r_list=[7, 10], ax=ax, color=BLUE,
        palette=['#6796C3', BLUE], r_text_ha='center'
    )

    # mesh r vs f
    ax = axs['h']
    p_list = np.array(ny_data.iloc[0]["p"])
    p = p_list[np.abs(p_list - 5e-4).argmin()]
    r_vs_f_lineplot(data=ny_data.query('r<=20'), p=p, ax=ax, color=BLUE)

    # mesh n vs f
    ax = axs['i']
    n_vs_f_lineplot(
        new_york_asymptotic_file_name, ax=ax, p=5e-4
    )

    # ---- third row ----
    # p
    ax = axs["j"]
    img = mpimg.imread(r"figures/comp_risk.png")
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(r"Component risk  $p$", fontsize=TEXT_SIZE + 2, color=BLUE)

    # R
    ax = axs['k']
    r_list = [1, 5, 10]
    graphs = list(ny_data.query('r in @r_list').iloc[:3].graph)
    plot_examples(graphs, ax, height_ratio=0.8)
    ax.set_title(r"Redundancy $R$", fontsize=TEXT_SIZE + 2, color=BLUE)

    # N
    asym_df = rw.read_nxjson(new_york_asymptotic_file_name).drop_duplicates('n')
    ax = axs['l']
    graphs = list(asym_df.sort_values('n').graph.iloc[[3, 5, 7]])    
    plot_examples(graphs, ax, height_ratio=0.8)
    ax.set_title(r"Scale $N$", fontsize=TEXT_SIZE + 2, color=BLUE)

    # save
    if save:
        plt.savefig(
            r'outputs\ny_plot\ny_plot.svg',
            bbox_inches='tight',
            transparent=True,
        )
    return fig, axs
        

def ring_p_vs_f_zoom(
    new_york_file_name: str = NEW_YORK_FILE_NAME,
    xmax: float = 1.5e-5,
    save: bool = False,
) -> tuple[plt.Figure, plt.Axes]:
    """Separate figure: ring network P vs F, simple lineplot zoomed to xmax on both axes."""
    ny_data = rw.read_nxjson(new_york_file_name)

    # prepare data (same logic as p_vs_f_lineplot)
    data_for_plot = explode_data(ny_data)
    data_for_plot = data_for_plot.query('r == 1').copy()
    data_for_plot = data_for_plot.query('p != 0').copy()
    data_for_plot.loc[len(data_for_plot)] = pd.Series({'p': 0, 'SAIDI': 0, 'r': 1})
    data_for_plot.sort_values('p', inplace=True)
    fig_zoom, ax = plt.subplots(figsize=(2 / cm_to_inch, 2 / cm_to_inch))


    # grey y = x reference line
    sns.lineplot(x=[0, xmax], y=[0, xmax], linestyle='--', color=GREY, ax=ax, linewidth=LINEWIDTH)

    # data line
    sns.lineplot(
        data=data_for_plot,
        x='p', y='SAIDI',
        color=BLUE,
        marker='o',
        linewidth=LINEWIDTH,
        ax=ax,
        legend=False,
    )

    # axes limits (square)
    axis_scale = 1.05
    ax.set_xlim(-xmax / 100, xmax * axis_scale)
    ax.set_ylim(-xmax / 100, xmax * axis_scale)

    # inset-style axis formatting
    INSET_LABEL_SIZE = LABEL_SIZE
    INSET_TICK_SIZE  = TICK_SIZE  
    SPINE_LW = LINEWIDTH

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(1.5)
    configure_axes(ax)

    # fixed tick labels: 0 and 1e-5 only
    ax.set_xticks([0, 1e-5])
    ax.set_yticks([0, 1e-5])
    ax.set_xticklabels(["0", r"$10^{-5}$"])
    ax.set_yticklabels(["0", r"$10^{-5}$"])
    ax.tick_params(axis='both', which='major', labelsize=INSET_TICK_SIZE, width=1.5, length=6)

    # pc point and line
    pc = find_pc(data_for_plot)
    if pc is not None and pc <= xmax:
        ax.axvline(x=pc, color=YELLOW, linestyle='--', linewidth=SPINE_LW)
        ax.scatter([pc], [pc], color=YELLOW, s=POINTSIZE, zorder=10)

    ax.set_xlabel("", fontsize=INSET_LABEL_SIZE)
    ax.set_ylabel("", rotation=90, fontsize=INSET_LABEL_SIZE)

    if save:
        fig_zoom.savefig(
            r'outputs\ny_plot\ny_plot_ring_zoom.svg',
            bbox_inches='tight',
            transparent=True,
        )
    return fig_zoom, ax
        

def make_figure_layout(
    figsize: tuple = (15, 15),
    big_ratio: float = 1.3,
    group_gap: float = 0.35,
    inter_group_gap: float = 0.3,
    row2_wspace: float = 0.35,
    row3_wspace: float = 0.35,
    top_hspace: float = 0.25,
    label_style: dict = None,
) -> tuple[plt.Figure, dict]:
    """
    Layout:
      Row 1: two groups, each group is 2x2 with left axis spanning two rows.
             big axes (a,d) have width big_ratio relative to small axes.
      Row 2: 3 equal columns (g,h,i)
      Row 3: 3 equal columns (j,k,l)

    Returns
    -------
    fig : matplotlib.figure.Figure
    ax : dict mapping letter -> Axes (e.g., ax['a'], ..., ax['l'])
    """
    if label_style is None:
        label_style = dict(fontsize=12, fontweight="bold")

    fig = plt.figure(figsize=figsize, constrained_layout=False)

    # Top-level: 3 equal-height rows
    gs_outer = fig.add_gridspec(
        nrows=3, ncols=1,
        height_ratios=[2, 1, 1],
        hspace=top_hspace
    )

    axes = {}

    # -------------------------
    # Row 1: nested 2x5 grid
    # cols: [bigL, smallL, gap, bigR, smallR]
    # -------------------------
    gs_top = gs_outer[0].subgridspec(
        nrows=2, ncols=7,
        width_ratios=[big_ratio, inter_group_gap, 1.0, group_gap, big_ratio, inter_group_gap, 1.0],
        wspace=0,
        hspace=0.5
    )

    # Left group
    axes["a"] = fig.add_subplot(gs_top[:, 0])   # spans 2 rows
    axes["b"] = fig.add_subplot(gs_top[0, 2])
    axes["c"] = fig.add_subplot(gs_top[1, 2])

    # Right group
    axes["d"] = fig.add_subplot(gs_top[:, 4])   # spans 2 rows
    axes["e"] = fig.add_subplot(gs_top[0, 6])
    axes["f"] = fig.add_subplot(gs_top[1, 6])

    # -------------------------
    # Row 2: 1x3 equal
    # -------------------------
    gs_mid = gs_outer[1].subgridspec(
        nrows=1, ncols=3,
        wspace=row2_wspace
    )
    axes["g"] = fig.add_subplot(gs_mid[0, 0])
    axes["h"] = fig.add_subplot(gs_mid[0, 1])
    axes["i"] = fig.add_subplot(gs_mid[0, 2])

    # -------------------------
    # Row 3: 1x3 equal
    # -------------------------
    gs_bot = gs_outer[2].subgridspec(
        nrows=1, ncols=3,
        wspace=row3_wspace
    )
    axes["j"] = fig.add_subplot(gs_bot[0, 0])
    axes["k"] = fig.add_subplot(gs_bot[0, 1])
    axes["l"] = fig.add_subplot(gs_bot[0, 2])

    # Label each axis: a, b, c, ...
    for letter, ax in axes.items():
        ax.text(
            -0.05, 1, f'{letter}', color='black', fontsize=TEXT_SIZE + 2, transform=ax.transAxes, ha='right'
        )

    return fig, axes


##### line graphs ####
def p_vs_f_lineplot(
    data: pd.DataFrame, r_list: list, ax: plt.Axes, color: str = GREEN, 
    r_text_ha: str = 'center', xmax: float = None, palette: list = None,
) -> None:
    """Plot P vs F lines for specified redundancy values."""
    # if per-r palette not given, fall back to the single color repeated
    _palette = palette if palette is not None else [color]
    # transform data
    data_for_plot = explode_data(data)
    pmax = data_for_plot["p"].max()
    xlim = pmax if xmax is None else xmax
    ylim = pmax
    # filter
    data_for_plot  = data_for_plot.query('r in @r_list').copy()
    # y=x line
    sns.lineplot(
        x=[0, xlim], y=[0, xlim], linestyle='--', color=GREY, ax=ax, linewidth=LINEWIDTH
    )
    # the saidi line
    data_for_plot = data_for_plot.query('p != 0').copy()
    for r in data_for_plot['r'].unique():
        data_for_plot.loc[len(data_for_plot)] = pd.Series({'p': 0, 'SAIDI': 0, 'r': r})
    data_for_plot.sort_values(['r', 'p'], inplace=True)
    sns.lineplot(
        data=data_for_plot,
        x='p',
        y='SAIDI',
        hue='r',
        palette=_palette,
        marker='o',
        linewidth=3,
        ax=ax
    )
    # add pc
    for r, group in data_for_plot.groupby('r'):
        if r <= 1:
            continue
        pc = find_pc(group.query('p > 1e-5'))
        print(f"____ pc = {pc} ____")
        if pc:
            ax.scatter(
                x=[pc], y=[pc], color=[YELLOW], s=POINTSIZE, zorder=10, alpha=1
            )
            ax.axvline(x=pc, color=YELLOW, linestyle='--', linewidth=LINEWIDTH)
            ax.text(
                x=pc - 0.00001, y=5e-4, s=r'$p_c$',
                va='bottom', ha='right', fontsize=TEXT_SIZE, color=YELLOW, zorder=100
            )
    # add r
    for i, (r, group) in enumerate(data_for_plot.query('r in @r_list').groupby('r')):
        max_p = group[group['SAIDI'] <= pmax]['p'].max()
        text_color = _palette[i % len(_palette)]
        ax.text(
            max_p, pmax, s=rf"$R={int(r)}$", ha=r_text_ha, va='bottom',
            fontsize=TEXT_SIZE, color=text_color,
            bbox=dict(facecolor="white", edgecolor="white", boxstyle="round,pad=0.3", alpha=1)
        )
    
    # add F=p
    ax.text(
        0.8, 0.8, r"$F=p$", rotation=0, color=GREY,
        ha='left', va='top', fontsize=TEXT_SIZE, transform=ax.transAxes
    )

    # config axes
    axis_scale = 1.05
    ax.set_xlim((-xlim / 100, xlim * axis_scale))
    ax.set_ylim((-xlim / 100, ylim * axis_scale))
    ax.legend().set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    configure_axes(ax, )
    set_custom_scientific_format(ax, axis='x', factor=int(np.floor(np.log10(xlim))), int_labels=False)
    set_custom_scientific_format(ax, axis='y', factor=int(np.floor(np.log10(ylim))), int_labels=False, )
    ax.set_xlabel(r"$p$", fontsize=LABEL_SIZE)
    ax.set_ylabel(r"$F$", rotation=90, fontsize=LABEL_SIZE)

def r_vs_f_lineplot(
    data: pd.DataFrame, p: float, ax: plt.Axes, color: str = GREEN
) -> None:
    """Plot Redundancy (R) vs F for a specific p value."""
    data_for_plot = explode_data(data)
    data_for_plot = data_for_plot.query("p == @p")
    # the saidi line
    sns.lineplot(
        data=data_for_plot,
        x='r',
        y='SAIDI',
        palette=[color],
        marker='o',
        linewidth=3,
        ax=ax
    )
    # add F=p line
    ax.axhline(
        y=p, color=RED, linestyle='--', linewidth=LINEWIDTH
    )
    # add rc
    rc = find_rc(data_for_plot, p)
    print(f"____ rc = {rc} ____")
    ax.scatter(
        x=[rc], y=[p], color=[YELLOW], s=POINTSIZE, zorder=10, alpha=1
    )
    ax.text(x=rc - 1, y=p - 0.0001, s=r'$R_c$', va='top', ha='center', fontsize=TEXT_SIZE, color=YELLOW)
    # add p label
    ax.text(
        data_for_plot['r'].max(), p, s=fr'$p=${format_p_scientific(p)}',
        ha='right', va='bottom', fontsize=TEXT_SIZE, color=RED
    )
    # configure axes
    ax.set_yscale('log')
    ax.legend().set_visible(False)
    configure_axes(ax)
    ax.set_xlabel(r"$R$", fontsize=LABEL_SIZE)
    ax.set_ylabel(r"$F$", rotation=90, fontsize=LABEL_SIZE)

def n_vs_f_lineplot(new_york_asymptotic_file_name: str, p: float, ax: plt.Axes) -> None:
    """Plot Network Size (N) vs F for a specific p value."""
    data = rw.read_nxjson(new_york_asymptotic_file_name)
    data_for_plot = explode_data(data)
    # line plot
    sns.lineplot(
        data=data_for_plot,
        x='n',
        marker='o',
        y='SAIDI',
        linewidth=3,
        ax=ax
    )
    # hline and text
    ax.axhline(p, color=RED, linestyle='--', linewidth=LINEWIDTH)
    max_n = data_for_plot['n'].max()
    ax.text(max_n * 0.95, p * 1.02, rf"$p=${format_p_scientific(p)}", ha='right', va='bottom', color=RED, fontsize=TEXT_SIZE)
    # plot Nc
    nc = find_nc(data_for_plot, p)
    print(f"____ nc = {nc} ______")
    if nc:
        ax.scatter(
            x=[nc], y=[p], color=[YELLOW], s=POINTSIZE, zorder=10, alpha=1
        )
        ax.text(x=nc * 0.9, y=p * 1.05, s=r'$N_c$', va='bottom', ha='right', fontsize=TEXT_SIZE, color=YELLOW)

    # configure axes
    ax.set_yscale('log')
    ax.set_xscale('log')
    configure_axes(ax)
    ax.set_xlabel(r'$N$', fontsize=LABEL_SIZE)
    ax.set_ylabel(r'$F$', rotation=90, fontsize=LABEL_SIZE)




##### helper #####
def explode_data(data: pd.DataFrame) -> pd.DataFrame:
    exploded = data.explode(["p", "rel"])
    if 'trial' not in exploded.columns:
        exploded['trial'] = 0
    exploded = exploded.groupby(["n", "m", "r", "p"])[["rel"]].apply('mean').reset_index()
    exploded.sort_values(['r', 'p'], inplace=True)
    exploded.rename(columns={'rel': 'SAIDI'}, inplace=True)
    exploded.reset_index(drop=True, inplace=True)
    return exploded

def max_val_below(x: np.array, thresh: float) -> float:
    x_filterd = x[x <= thresh]
    if len(x_filterd) == 0:
        x.sort()
        return x[1]
    return thresh

def find_pc(group: pd.DataFrame) -> float | None:
    """
    Find the first point where SAIDI = p (where p is from the data).
    Sorts by p and finds the intersection using linear interpolation.
    """
    group = group.sort_values('p').reset_index(drop=True).iloc[1:]
    p_vals = group['p'].values.astype(float)
    saidi_vals = group['SAIDI'].values.astype(float)
    
    # Find where SAIDI crosses p
    diff = saidi_vals - p_vals
    sign_changes = np.where(np.diff(np.sign(diff)) != 0)[0]
    
    if len(sign_changes) == 0:
        return None
    
    idx = sign_changes[0]
    # Direct linear interpolation: find x where y = 0
    # Using: x = x1 - y1 * (x2 - x1) / (y2 - y1)
    x1, x2 = p_vals[idx], p_vals[idx + 1]
    y1, y2 = diff[idx], diff[idx + 1]
    
    if abs(y2 - y1) < 1e-15:
        return (x1 + x2) / 2
    
    pc = x1 - y1 * (x2 - x1) / (y2 - y1)
    return pc

def find_rc(group: pd.DataFrame, p: float) -> float:
    """
    Find the first point where SAIDI = p (where p is a parameter).
    Sorts by r and finds the intersection at the given p value.
    """
    group = group.sort_values('r').reset_index(drop=True)
    r_vals = group['r'].values.astype(float)
    saidi_vals = group['SAIDI'].values.astype(float)
    
    # Find where SAIDI crosses the horizontal line at y=p
    diff = saidi_vals - p
    sign_changes = np.where(np.diff(np.sign(diff)) != 0)[0]
    
    if len(sign_changes) == 0:
        return None
    
    idx = sign_changes[0]
    # Direct linear interpolation: find x where SAIDI = p
    # Using: x = x1 - y1 * (x2 - x1) / (y2 - y1) where y = SAIDI - p
    x1, x2 = r_vals[idx], r_vals[idx + 1]
    y1, y2 = diff[idx], diff[idx + 1]
    
    if abs(y2 - y1) < 1e-15:
        return (x1 + x2) / 2
    
    rc = x1 - y1 * (x2 - x1) / (y2 - y1)
    return rc

def find_nc(group: pd.DataFrame, p: float) -> float:
    """
    Find the first point where F (SAIDI) = p (where p is a parameter).
    Sorts by n and finds the intersection at the given p value.
    Works for both increasing and decreasing data.
    """
    group = group.sort_values('n').reset_index(drop=True)
    n_vals = group['n'].values.astype(float)
    saidi_vals = group['SAIDI'].values.astype(float)
    
    # Find where SAIDI crosses the horizontal line at y=p
    diff = saidi_vals - p
    sign_changes = np.where(np.diff(np.sign(diff)) != 0)[0]
    
    if len(sign_changes) == 0:
        return None
    
    idx = sign_changes[0]
    # Direct linear interpolation: find x where SAIDI = p
    # Using: x = x1 - y1 * (x2 - x1) / (y2 - y1) where y = SAIDI - p
    x1, x2 = n_vals[idx], n_vals[idx + 1]
    y1, y2 = diff[idx], diff[idx + 1]
    
    if abs(y2 - y1) < 1e-15:
        return (x1 + x2) / 2
    
    nc = x1 - y1 * (x2 - x1) / (y2 - y1)
    return nc

def _key(e: tuple) -> tuple:
    return tuple(sorted(e))

def _minimum_spanning_tree(G: nx.Graph) -> nx.Graph:
    for e in G.edges:
        weight = np.linalg.norm(
            np.array(G.nodes[e[0]]["pos"]) -
            np.array(G.nodes[e[1]]["pos"])
        )
        G.edges[e]["weight"] = weight
    return nx.minimum_spanning_tree(G, weight="weight")




###### lower row #####

def plot_examples(
        graphs: list[nx.Graph],
        ax: plt.Axes,
        height_ratio: float = 1
) -> None:
    """Plot example networks side by side."""
    plot_graphs_side_by_side(graphs, ax)
    lower, upper = ax.get_ylim()
    r = (2 - height_ratio)
    new_upper = (1 - r) * lower + r * upper
    new_lower = r * lower + (1 - r) * upper
    ax.set_ylim(new_lower, new_upper)

    

def draw_convex_hull_background(
    ax,
    pos,
    *,
    pad=0.05,
):
    """
    Compute convex hull of points and draw it as a background polygon.
    """
    # ---- extract coordinates ----
    if isinstance(pos, dict):
        pts = np.asarray(list(pos.values()), dtype=float)
    else:
        pts = np.asarray(pos, dtype=float)

    # remove NaNs
    if pts.size == 0:
        return None
    pts = pts[np.isfinite(pts).all(axis=1)]
    if pts.shape[0] < 3:
        return None

    # ---- compute hull ----
    poly_xy = None

    # Try Shapely first (nice because we can buffer outward easily)
    from shapely.geometry import MultiPoint
    hull = MultiPoint(pts).convex_hull
    if hull.is_empty or hull.geom_type != "Polygon":
        return None
    if pad and pad != 0:
        hull = hull.buffer(pad)
        if hull.is_empty or hull.geom_type != "Polygon":
            return None
    poly_xy = np.asarray(hull.exterior.coords, dtype=float)

    from matplotlib.patches import Polygon
    ax.add_patch(Polygon(
        poly_xy, closed=True,
        facecolor='none',
        alpha=0.6,
        edgecolor='black',
        linewidth=1.5,
        clip_on=False,
        zorder=0
    ))
    return poly_xy

def net_plot_polygon(graph: nx.Graph, ax: plt.Axes, pad: float = 0.05) -> None:
    """Draw network on a polygon."""
    # draw network
    pos = nx.get_node_attributes(graph, "pos")
    draw_network(
        graph,
        pos=pos,
        node_color=GREY,
        node_size=0,
        edge_color="grey",
        with_labels=False,
        ax=ax
    )

    # Redundant edges
    T_edges = set(map(_key, _minimum_spanning_tree(graph).edges))
    re_edges = [e for e in graph.edges if _key(e) not in T_edges]
    draw_edges(
        graph,
        pos=pos,
        edgelist=re_edges,
        width=1,
        edge_color=[BLUE] * len(re_edges),
        ax=ax,
        zorder=100,
    )

def plot_graphs_side_by_side(
        graphs: list[nx.Graph],
        ax: plt.Axes
) -> None:
    """Plot multiple graphs side by side within the same axes."""
    all_points = np.array([data["pos"] for graph in graphs for _, data in graph.nodes.items()])
    xmin = all_points[:, 0].min()
    xmax = all_points[:, 0].max()
    ymin = all_points[:, 1].min()
    ymax = all_points[:, 1].max()

    k = len(graphs)
    def get_pos(pos: dict, i: int) -> dict:
        return {
            node: (
                (p[0] - xmin) / (xmax - xmin) / k + i / k,
                (p[1] - ymin) / (ymax - ymin)
            )
            for node, p in pos.items()
        }

    for i, graph in enumerate(graphs):
        # rescale
        pos = nx.get_node_attributes(graph, "pos")
        new_pos = get_pos(pos, i)
        new_graph = graph.copy()
        nx.set_node_attributes(new_graph, new_pos, "pos")
        net_plot_polygon(new_graph, ax, pad=0)
    # border
    last_graph = graphs[-1]
    pos = nx.get_node_attributes(last_graph, "pos")
    for i in range(k):
        draw_convex_hull_background(ax, get_pos(pos, i), pad=0.01)

    ax.set_xlim((0.0, 1.0))
    ax.set_ylim((0.0, 1.0))

##### MCMC ####
def reliability_MCMC_plot(
    r: float,
    data_file_name: str,
    ax: plt.Axes = None,
) -> None:
    """Plot reliability SAIDI metric over time from MCMC simulation."""
    # read data
    data = rw.read_nxjson(data_file_name)
    data = data[data["r"] == r].copy()
    if len(data) == 0:
        raise ValueError(f"data has nor r valuf of {r}")
    data.groupby("p").apply(resample_cum_saidi, dt=0.5).reset_index(level=0)

    # transform data
    data['idx'] = data.groupby('p')['t'].rank()
    data = data.sort_values(['p', 'idx'])
    data['p_str'] = data['p'].astype(str) 
    t_max = data['t'].max()
    colors = list(sns.color_palette('rocket_r'))
    
    # plot saidi
    ymin = data.query('t == @t_max')['saidi'].min() * 1e-2
    sns.lineplot(
        data=data.query('saidi > @ymin'),
        x='t',
        y='saidi',
        hue='p',
        palette=colors,
        ax=ax,
        linewidth=3,
        zorder=10,
        legend=False
    )
   
    # add mean saidi line
    for i, (p, group) in enumerate(data.groupby('p')):
        saidi = group.sort_values('t').iloc[-1]['saidi']
        ax.plot(
            [0, t_max], [saidi, saidi], linestyle='--', color=colors[i], linewidth=2
        )
        ax.text(
            t_max * 0.4, saidi * 1.05, r"$p=$" + str(p),
            color=colors[i], ha='left', va='bottom', fontsize=ANNOTATION_SIZE
        )
        i += 1
    
    # configure
    ax.set_xlim((5, t_max* 1.3))
    ax.set_yscale('log')
    ax.set_xticks(list(range(0, int(t_max) + 1, 100))[1:])
    ax.set_ylim([ymin, ax.get_ylim()[1]])
    ax.set_xlabel('days', fontsize=LABEL_SIZE)
    ax.set_ylabel(r'$F$', rotation=90, fontsize=LABEL_SIZE)
    configure_axes(ax)


def resample_cum_saidi(df: pd.DataFrame, dt: float):
    """
    Resample cumulative SAIDI to constant dt using linear interpolation.
    """
    p = df.iloc[0]["p"]
    df.loc[len(df)] = {'t': 0, 'saidi': 0}
    df.sort_values('t', inplace=True)
    t_end = df["t"].iloc[-1]

    t_uniform = np.arange(0.0, t_end + dt, dt)

    saidi_uniform = np.interp(
        t_uniform,
        df["t"].to_numpy(),
        df["saidi"].to_numpy(),
    )

    return pd.DataFrame({
        "t": t_uniform,
        "saidi": saidi_uniform,
    })

def style_sns_legend(
    ax,
    *,
    title_fontsize=TEXT_SIZE,
    text_fontsize=TEXT_SIZE,
    reverse=False
):
    leg = ax.get_legend()
    if leg is None:
        return

    # --- extract current legend entries ---
    handles = leg.legend_handles
    labels = [t.get_text() for t in leg.get_texts()]
    title = leg.get_title().get_text()

    # remove old legend
    leg.remove()

    # reverse order if requested
    if reverse:
        handles = handles[::-1]
        labels = labels[::-1]

    # --- rebuild legend ---
    leg = ax.legend(
        handles,
        labels,
        title=title,
        loc="lower right",
        bbox_to_anchor=(1.0, 0),
        frameon=False,
    )

    # --- styling ---
    leg.get_title().set_fontsize(title_fontsize)
    leg.get_title().set_ha("left")

    for txt in leg.get_texts():
        txt.set_fontsize(text_fontsize)

    # left-align whole legend block (important!)
    try:
        leg._legend_box.align = "left"
    except Exception:
        pass
