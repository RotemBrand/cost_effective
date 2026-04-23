import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from indexes.probs import Poly, Array
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from indexes.utilities import create_sparse_graph
import networkx as nx
import figures.ny_simulation as NY
import utilities.read_write as rw
from tqdm import tqdm
from pathlib import Path
from figures.ny_plot import NEW_YORK_ASYMPTOTIC_FILE_NAME
from utilities.figures_utilities import (
    RED, GREEN, GREY, BLUE, YELLOW,
    TEXT_SIZE, LABEL_SIZE, TICK_SIZE, ANNOTATION_SIZE,
    LINE_WIDTH, cm_to_inch, configure_axes, format_p_scientific,
    get_xylim_from_graphs, set_custom_scientific_format
)
from indexes import GraphRel
from utilities import draw_network

tqdm.pandas()

DATA_FOLDER = Path("data/topologies")
MARKER_WIDTH = 1.5





############ plot ############
def nature_topologies_graph(
        trees_file: Path = DATA_FOLDER / "trees.nxjson",
        rings_file: Path = DATA_FOLDER / "rings.nxjson",
        rings_poly_file: Path = DATA_FOLDER / "poly_data.nxjson",
        regular_file: Path = DATA_FOLDER / "regulars.nxjson",
        save: bool = False,
) -> np.ndarray:
    """
    Generate the main topologies figure.
    
    Creates a 3x3 subplot grid analyzing tree, ring, and 3-regular graph topologies.
    """
    # read data
    trees = rw.read_nxjson(trees_file)
    rings = rw.read_nxjson(rings_file)
    rings_poly = rw.read_nxjson(rings_poly_file)
    print(rings_poly.index)
    regulars = rw.read_nxjson(regular_file)
    # create figure
    fig, axs = plt.subplots(
        3, 3, figsize=(3*10/ cm_to_inch, 3*10/ cm_to_inch), gridspec_kw={'hspace': 0.4, 'wspace': 0.5}
    )

    # trees plots
    plot_idx = [87, 9, 77]
    trees_example(ax=axs[0, 0])
    plot_trees_rel(trees, ax=axs[0, 1], plot_idx=plot_idx)
    plot_avg_d_rel(trees, ax=axs[0, 2])

    # rings plots

    all_n = np.array(sorted(set(rings_poly.n)))
    n_list = all_n[[-1, -2, -3]]
    rings_example(axs[1, 0])
    plot_ring_rel(rings_poly, n_list=n_list, pmax=5e-3, ax=axs[1, 1])
    plot_pc_by_n(rings_poly, ax=axs[1, 2])

    # regulars plots
    regulars_example(axs[2, 0])
    all_n = np.array(sorted(set(regulars['n'])))
    n_list = all_n[[-5, -3, -1]]
    plot_regular_rel(regulars=regulars, ax=axs[2, 2], n_list=n_list)
    plot_regular_rel_by_n(regulars=regulars, ax=axs[2, 1])

    # numbers
    number_plots(axs)

    # save
    if save:
        plt.savefig(
            r'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\topologies_plot\topologies_plot.svg',
            bbox_inches='tight',
            transparent=True,
            format='svg',
        )
    # trees
    save_selected_trees(trees, plot_idx=plot_idx, save=save)
    return axs

def empty_axes_with_color(ax: plt.Axes, color: str, title: str) -> None:
    """Clear an axis and set a colored title on the y-axis."""
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel(title, fontsize=TEXT_SIZE, color=color)
    ax.set_frame_on(True)
    ax.yaxis.set_label_position('left')

LETTERS = ['a', 'b', 'c' ,'d', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r']

def number_plots(axs: np.ndarray, start: int = 0) -> None:
    """Add letter labels (a, b, c...) to each subplot."""
    for i, ax in enumerate(axs.flatten()):
        ax.text(
            -0.04, 1.0, f"{LETTERS[start + i]}",
            fontsize=TEXT_SIZE, color="black",
            ha='right', va='bottom',
            transform=ax.transAxes,
        )


###### saidi ####
def get_saidi(graph: nx.Graph, sources: list[int], p_list: list[float], rng: np.random.Generator, mean_cycle_days: float = 0.5) -> list[float]:
    """Calculate SAIDI values for a graph over a list of failure probabilities."""
    edges_prob = {e: Array(p_list) for e in graph.edges}
    gr = GraphRel(graph, sources=sources, edges_prob=edges_prob)
    res = gr.calc_rel_simulation(T_days=365*5, mean_cycle_days=mean_cycle_days, rng=rng, show_progress=False)
    saidi = [r.rel_result for r in res]
    return saidi

###### trees #####
### plot ###
def plot_trees_rel(trees: pd.DataFrame, plot_idx: list, ax: plt.Axes = None) -> None:
    """Plot reliability F vs p for various tree topologies."""
    if ax is None:
        ax = plt.gca()
    # plot trees rel points
    upper = 1e-2 * 1.2
    trees['size'] = pd.Series(trees.index.isin(plot_idx)).map({True: 150, False: 50})
    trees_ = trees.query('saidi < @upper').sort_values('size')
    sns.scatterplot(
        data=trees_,
        x='p',
        y='saidi',
        color=RED,
        s=trees_['size'],
        marker='o',
        hue_order='size',
        facecolor=['white' if idx in plot_idx else 'none' for idx in trees_.index],
        edgecolor=RED,
        linewidth=MARKER_WIDTH,
        ax=ax
    )
    sns.lineplot(
        x=[0, upper],
        y=[0, upper],
        color=GREY,
        linestyle='--',
        linewidth=LINE_WIDTH - 1,
        ax=ax
    )
    # configure axes
    ax.legend().set_visible(False)
    ax.set_xlim((0, upper))
    ax.set_ylim((0, upper))
    ax.set_ylabel(r"$F$")
    ax.set_xlabel(r"$p$")
    set_custom_scientific_format(ax, axis='both', factor=-3)
    configure_axes(ax)

def save_selected_trees(trees: pd.DataFrame, plot_idx: list, save: bool = False) -> None:
    """Save selected tree graphs as SVG files."""
    selected_trees = trees.loc[plot_idx]
    selected_trees.sort_values('p', inplace=True)
    for i, (_, row) in enumerate(selected_trees.iterrows()):
        G = row.tree
        fig, ax = plt.subplots(figsize=(2.5/ cm_to_inch, 2.5/ cm_to_inch))
        draw_network(
            G,
            ax=ax,
            node_color=[RED if node == row.source else 'white' for node in G.nodes],
            with_labels=False,
            edge_color=RED,
            edgecolors=RED,
            node_size=30,
        )
        ax.axis(True)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linewidth(1)
            spine.set_capstyle('round')
            spine.set_joinstyle('round')
        if save:
            fig.savefig(
                fr'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\topologies_plot\tree_{i}.svg',
                bbox_inches='tight',
                transparent=True,
                format='svg',
            )

        

def plot_avg_d_rel(trees: pd.DataFrame, ax: plt.Axes = None) -> None:
    """Plot average distance * p vs SAIDI (F)."""
    if ax is None:
        ax = plt.gca()
    # plot avg_d*p vs saidi
    trees['avg_d*p'] = trees['avg_d'] * trees['p']
    sns.scatterplot(
        data=trees,
        x='avg_d*p',
        y='saidi', 
        color=RED,
        facecolor='none',
        edgecolor=RED,
        linewidth=MARKER_WIDTH,
        ax=ax
    )
    # plot y=x line
    lower = trees[['saidi', 'avg_d*p']].min().min()
    upper = trees[['saidi', 'avg_d*p']].max().max()
    sns.lineplot(
        x=[lower, upper],
        y=[lower, upper],
        color=GREY,
        linestyle='--',
        linewidth=LINE_WIDTH - 1,
        ax=ax
    )
    ax.set_ylabel(r"$F$")
    ax.set_xlabel(r"$F_{Theory}$")
    set_custom_scientific_format(ax, axis='both', factor=-2)
    configure_axes(ax)

def trees_example(ax: plt.Axes) -> None:
    """Plot example tree graphs as inset axes."""
    # generate trees
    star = nx.star_graph(8)
    path = nx.path_graph(6)
    binary_tree = nx.balanced_tree(r=2, h=3)
    np.random.seed(42)
    prufer = np.random.randint(0, 16, 14)
    generic_tree = nx.from_prufer_sequence(prufer)
    # configure plot
    empty_axes_with_color(ax, color=RED, title=r"$G_{Tree}$")
    graphs = [star, binary_tree, path, generic_tree]
    inset_grid = []
    # create the inset grid
    for i in range(2):
        for j in range(2):
            inset_ax = inset_axes(ax, width="90%", height="90%", loc='center',
                                  bbox_to_anchor=(j * 0.5, 1 - (i + 1) * 0.5, 0.5, 0.5),
                                  bbox_transform=ax.transAxes, borderpad=0)
            inset_grid.append(inset_ax)
    # plot graphs
    i = 0
    for g, sub_ax in zip(graphs, inset_grid):
        draw_network(
            g,
            ax=sub_ax,
            node_color=[RED if node == 0 else 'white' for node in g.nodes],
            with_labels=False,
            edge_color=RED,
            edgecolors=RED,
            node_size=40,
            arc_edges=i == 0
        )
        i += 1

### data ###

def generate_tree_data(max_size: int, n_trees: int, seed: int = 42) -> pd.DataFrame:
    """Generate and save simulated reliability data for random trees."""
    trees = create_random_trees(max_size=max_size, n_trees=n_trees, seed=seed)
    trees = simulate_trees_rel(trees)
    # save
    rw.write_nxjson(trees, DATA_FOLDER / "trees.nxjson")
    return trees

def random_dist_from_source(tree: nx.Graph, source: int) -> float:
    """
    Returns the average shortest path length from the source to all nodes in the tree.
    """
    lengths = nx.single_source_shortest_path_length(tree, source)
    return np.mean(list(lengths.values()))

def create_random_trees(max_size: int, n_trees: int, seed: int = 42) -> pd.DataFrame:
    """Create a dataframe of random tree structures of various types."""
    data = []
    rng = np.random.default_rng(seed)
    for i in range(n_trees):
        size = rng.integers(5, max_size + 1)
        tree_type = i % 3
        if tree_type == 0:
            # Near-star: most entries are the same node, a few are random
            center = rng.integers(0, size)
            prufer = [center] * (size - 2)
            # Randomly replace a few entries
            for _ in range(max(1, (size - 2) // 5)):
                idx = rng.integers(0, size - 2)
                prufer[idx] = rng.integers(0, size)
        elif tree_type == 1:
            # Near-chain: random permutation, a few entries replaced with random nodes
            perm = rng.permutation(size)
            prufer = perm[:size - 2].tolist()
            for _ in range(max(1, (size - 2) // 5)):
                idx = rng.integers(0, size - 2)
                prufer[idx] = rng.integers(0, size)
        else:
            # Fully random
            prufer = rng.integers(0, size, size - 2).tolist()
        tree = nx.from_prufer_sequence(prufer)
        p = rng.uniform(1e-3, 1e-2)
        data.append((tree, p))
    data = pd.DataFrame(data, columns=['tree', 'p'])
    data['n'] = data['tree'].apply(lambda x: x.number_of_nodes())
    data['source'] = 0
    data['avg_d'] = data.apply(lambda row: random_dist_from_source(row['tree'], row['source']), axis=1)
    return data

def simulate_trees_rel(trees: pd.DataFrame) -> pd.DataFrame:
    """Simulate reliability SAIDI for a dataframe of trees."""
    rng = np.random.default_rng(356)
    def simulate_saidi(tree, p, source) -> float:
        return get_saidi(tree, sources=[source], p_list=[p], rng=rng)[0]
    trees['saidi'] = trees.apply(lambda row: simulate_saidi(*row[["tree", "p", "source"]]), axis=1)
    return trees


###### rings #####
### plot ###
def plot_ring_rel(poly_data: pd.DataFrame, n_list: list, ax: plt.Axes, pmax: float = 5e-3) -> None:
    """Plot reliability F vs p for ring topologies and their polynomial fits."""

    # plot the simulated points
    palette = sns.color_palette("Greens", len(n_list) + 1)[1:]
    print(n_list)
    poly_data_for_plot = poly_data.query('n in @n_list').drop_duplicates('n').set_index('n')
    # plot the fit of the real saidi
    fited_saidi = pd.DataFrame({'n': n_list})
    p_list = np.linspace(0, pmax, 100)
    fited_saidi['p'] = [p_list for _ in range(len(fited_saidi))]
    fited_saidi['saidi'] = fited_saidi.apply(
        lambda row: Poly(poly_data_for_plot.coeffs.loc[row.n])(row.p),
        axis=1
    )
    fited_saidi = fited_saidi.explode(['p', 'saidi'])
    fited_saidi = fited_saidi.query('saidi < @pmax')
    sns.lineplot(
        data=fited_saidi,
        x='p',
        y='saidi',
        hue='n',
        palette=palette,
        linewidth=LINE_WIDTH,
        linestyle='-',
        ax=ax
    )
    # plot the real ring saidi
    # plot F=p
    sns.lineplot(
        x=[0, pmax],
        y=[0, pmax],
        color='grey',
        linestyle='--',
        linewidth=LINE_WIDTH - 1,
        zorder=8,
        ax=ax
    )
    # plot pc
    for i, (n, row) in enumerate(poly_data_for_plot.iterrows()):
        sns.scatterplot(
            x=[row.pc],
            y=[row.pc],
            color=YELLOW,
            edgecolors=YELLOW,
            ax=ax,
            zorder=10
        )
        ax.axvline(row.pc, linestyle='--', color=YELLOW)
        if i == 2:
            ax.text(row.pc, row.pc, r"$p_c$", ha='right', va='bottom', fontsize=TEXT_SIZE, color=YELLOW)
    # plot labels
    # plot labels
    for i, n in enumerate(n_list[::-1]):
        ax.text(
            0.6, 0.1 + 0.1 * i, rf'$N={{{n}}}$',
            color=palette[i],
            ha='left', va='bottom',
            fontsize=TEXT_SIZE - 2,
            transform=ax.transAxes,
            bbox=dict(
                facecolor='white', edgecolor='none', boxstyle='round,pad=0', alpha=0.6
            )
        )
    ax.set_xlim((0, 3.8e-3))
    ax.set_ylim((0, 3.8e-3))
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$F$')
    ax.get_legend().set_visible(False)
    set_custom_scientific_format(ax, axis='both', factor=-3)
    configure_axes(ax)

def plot_pc_by_n(poly_data: pd.DataFrame, ax: plt.Axes = None) -> None:
    """Plot critical probability p_c versus network size N for rings."""
    if ax is None:
        ax = plt.gca()
    # plot pc points
    sns.scatterplot(
        data=poly_data,
        x='n',
        y='pc',
        color=GREEN,
        marker='o',
        facecolor='white',
        edgecolor=GREEN,
        linewidth=MARKER_WIDTH,
        ax=ax,
        legend=False,
        zorder=1000
    )
    # plot expected pc
    n = poly_data.n.values.astype('float')
    pc = poly_data.pc.values.astype('float')
    C = float(np.exp(np.mean(np.log(pc) + 2*np.log(n))))
    print(f"C={C:.4f}")
    expected = C * n ** (-2)
    sns.lineplot(
        x=n,
        y=expected,
        linestyle='-',
        color=GREEN,
        linewidth=LINE_WIDTH,
        ax=ax,
        label=r'$\sim N^{-2}$'
    )
    ax.legend(frameon=False)
    # configure axes
    ax.set_xlim((5, 2e2))
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$N$')
    ax.set_ylabel(r'$p_c$')

    configure_axes(ax)

def rings_example(ax: plt.Axes) -> None:
    """Plot example ring graphs as inset axes."""
    # create graphs to plot
    ring = nx.cycle_graph(10)
    strc = nx.complete_graph(4)
    conn2 = create_sparse_graph(strc, {e: sum(e) for e in strc.edges})
    empty_axes_with_color(ax, color=GREEN, title=r"$G_{Ring}$")
    # create two inset columns
    axes = []
    for j, graph in enumerate([conn2, ring]):
        inset_ax = inset_axes(ax, width="90%", height="100%", loc='center',
                             bbox_to_anchor=(j * 0.5, 0.25, 0.5, 0.5),
                             bbox_transform=ax.transAxes, borderpad=0)
        axes.append(inset_ax)
        draw_network(
            graph,
            ax=inset_ax,
            node_color=[GREEN if node == 0 else 'white' for node in graph.nodes],
            with_labels=False,
            edge_color=GREEN,
            edgecolors=GREEN,
            node_size=40,
        )

    
    
    

### data ###
def generate_rings_data(n_list: list, p_num: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate and save simulated reliability data for ring topologies."""
    # generate trees data
    rings = simulate_rings_criticality(n_list, p_num)
    # Fit degree 2 polynomial for each n
    poly_data = pd.DataFrame(columns=['coeffs'], index=rings['n'].unique())
    poly_data.index.name = 'n'
    for n, g in rings.groupby('n'):
        x = g['p'].values.astype(float)
        y = g['saidi'].values.astype(float)
        b2 = np.sum(y) / np.sum(x ** 2)
        coeffs = [0, 0, b2]
        poly_data.loc[n] = [coeffs]
    # Find intersection point for each n
    poly_data['pc'] = poly_data['coeffs'].apply(
        lambda coeff: get_y_is_x_point(Poly(coeff)) if coeff is not None
        else None
    )
    poly_data['n'] = poly_data.index
    # save
    rw.write_nxjson(rings, DATA_FOLDER / "rings.nxjson")
    rw.write_nxjson(poly_data, DATA_FOLDER / "poly_data.nxjson")
    return rings, poly_data

def simulate_rings_criticality(n_list: list, p_num: int) -> pd.DataFrame:
    """Simulate ring reliability for a range of probabilities to find criticality."""
    data = []
    rng = np.random.default_rng(1000)
    for n in tqdm(n_list):
        ring = nx.cycle_graph(n)
        p_min = float(n) ** (-2)
        p_max = 2 * float(n) ** (-1)
        p_list = np.logspace(np.log10(p_min), np.log10(p_max), p_num)
        saidi = get_saidi(ring, sources=[0], p_list=p_list, rng=rng)
        data.append([ring, n, p_list, saidi])
    data = pd.DataFrame(data, columns=['graph', 'n', 'p', 'saidi'])
    data_exploded = data.explode(['p', 'saidi'])
    return data_exploded

def get_y_is_x_point(coeff: list) -> float | None:
    """
    Find the point where the polynomial intersects the line y=x.
    """
    # Poly expects coefficients in increasing order (constant first)
    poly = Poly(coeff)
    p_line = Poly([0, 1])
    roots = (poly - p_line).roots()
    roots = roots[np.isreal(roots)].real
    roots = roots[roots > 0]
    if len(roots) == 0:
        return None
    return roots.max()


###### 3reg #####
### plot ###
def plot_regular_rel(regulars: pd.DataFrame, n_list: list, ax: plt.Axes = None) -> None:
    """Plot reliability F vs p for 3-regular graph topologies."""
    if ax is None:
        ax = plt.gca()

    # plot the simulated points
    palette = sns.color_palette("Blues", 4)[1:][::-1]
    data_for_plot = (
        regulars.query('n in @n_list')[["n", "p", "F"]]
        .explode(["p", "F"])
        .sort_values(["n", "p"]).iloc[::-1]
    )
    pmax = data_for_plot["p"].max()
    sns.lineplot(
        data=data_for_plot,
        x='p',
        y='F',
        hue='n',
        hue_order=n_list[::-1],
        marker='o',
        palette=palette,
        linewidth=LINE_WIDTH,
        linestyle='-',
        ax=ax
    )
    # plot F=p
    sns.lineplot(
        x=[0, pmax],
        y=[0, pmax],
        color='grey',
        linestyle='--',
        linewidth=LINE_WIDTH - 1,
        zorder=8,
        ax=ax
    )
    # plot labels
    for i, n in enumerate(n_list):
        ax.text(
            0.5, 0.23 + 0.1 * i, rf'$N={{{n}}}$',
            color=palette[::-1][i],
            ha='left', va='bottom',
            fontsize=TEXT_SIZE - 2,
            transform=ax.transAxes
        )
    # configure axes
    ax.set_xlim((0, pmax))
    ax.set_ylim((-pmax/10, pmax))
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$F$')
    ax.get_legend().set_visible(False)
    configure_axes(ax)

def plot_regular_rel_by_n(regulars: pd.DataFrame, ax: plt.Axes = None) -> None:
    """Plot reliability F versus network size N for 3-regular graphs."""
    if ax is None:
        ax = plt.gca()
    
    # plot the simulated points
    palette = sns.color_palette("Blues", 4)[1:]
    p = 5e-4
    data_for_plot = (
        regulars[["n", "p", "F"]]
        .explode(["p", "F"])
        .query('p == @p')
    )
    if len(data_for_plot) == 0:
        raise ValueError(f"p {p} is not a part of regulars df")
    sns.lineplot(
        data=data_for_plot,
        x='n',
        y='F',
        palette=palette,
        linewidth=LINE_WIDTH,
        linestyle='-',
        marker='o',
        ax=ax
    )
    # plot F=p
    p = 5e-4
    ax.axhline(p, color=RED, linestyle='--', linewidth=LINE_WIDTH - 1)
    # add p label
    ax.text(
        data_for_plot["n"].max(), p,
        s=fr'$p=${format_p_scientific(p)}',
        ha='right', va='bottom',
        fontsize=TEXT_SIZE, color=RED,
    )
    # configure axes
    ax.set_ylim((-1e-4, 1e-3))
    ax.set_xlabel(r'$N$')
    ax.set_ylabel(r'$F$')
    ax.set_xscale('log')
    set_custom_scientific_format(ax, axis='y', factor=-4)
    configure_axes(ax)

def regulars_example(ax: plt.Axes) -> None:
    """Plot example 3-regular graphs as inset axes."""
    # create graphs to plot
    regular = nx.random_regular_graph(d=3, n=30, seed=803)
    # plot graph
    inset_ax = inset_axes(ax, width="100%", height="100%", loc='center',
                         bbox_to_anchor=(0, 0, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
    draw_network(
        regular,
        ax=inset_ax,
        node_color=[BLUE if node == 20 else 'white' for node in regular.nodes],
        with_labels=False,
        edge_color=BLUE,
        edgecolors=BLUE,
        node_size=40,
        arc_edges=True
    )
    empty_axes_with_color(ax, color=BLUE, title=r"$G_{3Reg}$")

### data ###
def generate_random_regular_data(
        seed: int = 42, new_york_asymptotic_file_name: str = NEW_YORK_ASYMPTOTIC_FILE_NAME
    ) -> pd.DataFrame:
    """Generate and save simulated reliability data for 3-regular networks."""
    ny_data = rw.read_nxjson(new_york_asymptotic_file_name)
    ny_data.drop_duplicates("n", inplace=True)

    # improve
    def improve_to_3regular(graph: nx.Graph) -> nx.Graph:
        r = int(0.5 * len(graph) + 1)
        rng = np.random.default_rng(100)
        res = NY.random_improve_network(graph, r_list=[r], rng=rng)
        return res[r]
    ny_data['graph'] = [improve_to_3regular(graph) for graph in tqdm(ny_data["graph"], "improve")]

    # add rel
    rng = np.random.default_rng(2000)
    regulars = ny_data[["graph", "n", "sources"]].copy()
    p_list = list(np.linspace(0, 0.2, 10)) + [5e-4]
    polies = [
        get_saidi(row.graph, sources=row.sources, p_list=p_list, rng=rng, mean_cycle_days=10)
        for _, row in tqdm(list(regulars.iterrows()), "reliability")
    ]
    p = [list(p_list) for _ in polies]
    F = [list(poly) for poly in polies]
    regulars['p'] = p
    regulars['F'] = F

    # save
    rw.write_nxjson(regulars, DATA_FOLDER / "regulars.nxjson")
    return regulars



