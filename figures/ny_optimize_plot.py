import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import optimal_network as ON
from indexes.graph_rel import GraphRel
from figures.ny_plot import NEW_YORK_ASYMPTOTIC_FILE_NAME
import utilities.read_write as rw
import os
from tqdm import tqdm
from shapely.geometry import Point
from utilities.figures_utilities import (
    RED, GREEN, BLUE, GREY, cm_to_inch, configure_axes, format_p_scientific,
    LABEL_SIZE, TEXT_SIZE, TICK_SIZE, LINE_WIDTH,
    net_plot, saidi_with_lengths,
)

TEXT_SIZE = 17
LABEL_SIZE = 17
TICK_SIZE = 10
ANNOTATION_SIZE = 16
LINE_WIDTH = 2
MARKER_WIDTH = 1.5

# ===== plots =====
def ny_optimize_plot(nyc_optimal_path: str=r"data/nyc_optimal.nxjson", save: bool=False):
    fig, axs = plt.subplots(
        3, 2, figsize=(2 * 7 / cm_to_inch, 3 * 7 / cm_to_inch),
        gridspec_kw={'wspace': 0.4, 'hspace': 0.4}
    )
    # read data
    data = rw.read_nxjson(nyc_optimal_path)
    data = data.explode(["p", "f"]) # ["n", "alpha", "r", "p", "f", "graph"]
    pc_data = calculate_pc_data(data)

    # line plots
    lineplots(pc_data, axs)
    # draw netorks
    draw_sample_graph(data, n_idx=7, axs=axs)
    # number plots
    number_plots(axs)
    # save
    if save:
        plt.savefig(
            r'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\ny_optimize_plot\ny_optimize_plot.svg',
            bbox_inches='tight',
            dpi=300,
            transparent=True
        )
        plt.savefig(
            r'C:\Users\rotem\Desktop\מסמכים\תואר\תזה\write\baruch\ny_optimize_plot\ny_optimize_plot.pdf',
            bbox_inches='tight',
            dpi=300,
            transparent=True
        )
    return fig


def draw_sample_graph(data: pd.DataFrame, n_idx: int, axs):
    all_n = sorted(set(data['n']))
    all_alpha = sorted(set(data['alpha']))
    all_p = sorted(set(data['p']))
    n = all_n[n_idx]
    p = all_p[5]
    alpha_rows = {
        alpha: data.query('alpha == @alpha and n == @n and p == @p').iloc[0]
        for alpha in sorted(set(data['alpha']))
    }
    for alpha, ax in zip(all_alpha, axs[:, 0]):
        row = alpha_rows[alpha]
        graph = row.graph
        net_plot(graph, ax=ax, sources=row.sources, pad=[0.05, 0.05, 0.05, 0.3], basemap=True)
        ax.text(
            0.0, 0.95, fr'$\alpha={{{alpha:.1f}}}$', transform=ax.transAxes,
            ha='left', va='top', fontsize=TEXT_SIZE - 2, bbox=None
        )
        saidi = format_p_scientific(row.f)[1:-1]
        ax.text(
            0.0, 0.79, fr'$F={{{saidi}}}$', transform=ax.transAxes,
            ha='left', va='top', fontsize=TEXT_SIZE - 2
        )





def lineplots(pc_data: pd.DataFrame, axs):
    for (alpha, group), ax in zip(pc_data.groupby('alpha'), axs[:, 1]):
        sns.lineplot(
            data=group, x='n', y='pc_predicted', color='black', linewidth=LINE_WIDTH,ax=ax
        )
        sns.scatterplot(
            data=group, x='n', y='pc',
            edgecolor=BLUE, facecolor='white',
            linewidth=MARKER_WIDTH, ax=ax, zorder=100
        )
        beta = 2 * (1 - alpha)
        ax.text(
            0.9, 0.9, fr'$\beta={{{beta}}}$', ha='right', va='top',
            fontsize=TEXT_SIZE, transform=ax.transAxes
        )
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$N$', fontsize=LABEL_SIZE)
        ax.set_ylabel(r'$p_c$', fontsize=LABEL_SIZE)
        configure_axes(ax, label_size=LABEL_SIZE, tick_size=TICK_SIZE)

LETTERS = ['a', 'b', 'c' ,'d', 'e', 'f']
def number_plots(axs: np.ndarray):
    for i, ax in enumerate(axs.flatten()):
        ax.text(
            -0.03, 1.05, f"{LETTERS[i]}",
            fontsize=TEXT_SIZE, color="black",
            ha='right', va='bottom',
            transform=ax.transAxes,
        )



# ===== simulations =====
def simulate_ny_optimal_networks(
        nyc_asymptotic_path: str=NEW_YORK_ASYMPTOTIC_FILE_NAME,
        nyc_optimal_path: str=r"data/nyc_optimal.nxjson",
        skip_existing_data: bool=True
    ):
    # read data
    alpha_list = [0.3, 0.4, 0.5]
    p_list = np.logspace(-4, -2, 10)
    df = rw.read_nxjson(nyc_asymptotic_path).query('trial == 0')
    if skip_existing_data and os.path.exists(nyc_optimal_path):
        data = rw.read_nxjson(nyc_optimal_path)
    else:
        data = pd.DataFrame(columns=["r", "alpha", "n", "p", "f", "graph", "sources"])
    data.set_index(['alpha', 'n'], inplace=True)
    # construct optimal networks
    for alpha in alpha_list:
        i = 0
        rng = np.random.default_rng(100)
        for _, row in df.iterrows():
            n = row.n
            graph = row.graph
            r = int(n ** alpha)
            # check if data exist
            print(f'n = {n}, r = {r}, alpha = {alpha}')
            if (alpha, n) in data.index:
                continue
            points = np.array(list(nx.get_node_attributes(graph, "pos").values()))
            print("=== optimal network ===")
            strc_exact_vertices = (n > 1200) or (alpha == 0.5 and n > 500)
            graph = ON.optimal_network_from_points(
                points=points,
                r=r,
                kmeans_max_iter=7,
                seed=7,
                strc_n_init_iters=2,
                strc_exact_vertices=strc_exact_vertices,
                debug=False
            )
            print("=== simulate saidi ===")
            saidi_list = saidi_with_lengths(
                graph,
                sources=row.sources,
                p=p_list,
                mode="mean",
                T_days=365*5,
                mean_cycle_days=0.5,
                rng=rng,
                show_progress=True
            )
            data.loc[(alpha, n), :] = {
                'r': r, "p": p_list, 'f': saidi_list,
                'graph': graph, "sources": row.sources
            }
            print(f"{nyc_optimal_path=}")
            rw.write_nxjson(data.reset_index(), nyc_optimal_path)
            i += 1
    return data


def key(e):
    return tuple(sorted(e))

def read_sym_data(nyc_optimal_path: str) -> pd.DataFrame:
    data_json = rw.dict_read_json(nyc_optimal_path)
    data_json
    data = pd.DataFrame([
        [n, alpha, n_data["r"], p, f, n_data['graph']] 
        for alpha, a_data in data_json.items()
        for n, n_data in a_data.items()
        for p, f in n_data["saidi"].items()
    ], columns=["n", "alpha", "r", "p", "f", "graph"])
    return data

def calculate_pc_data(data: pd.DataFrame) -> pd.DataFrame:
    pc_data = []
    def get_pc(group):
        group = group.sort_values('p')
        group = group.iloc[5:]
        c = np.mean(group['f'] / (group['p'] ** 2))
        p_c = 1 / c
        return p_c
    # def get_pc_predicted(alpha_group: pd.DataFrame) -> float:
    #     alpha_group['pc_predicted'] = alpha_group['n'] ** (-2 * (1- alpha_group['alpha']))
    #     c = (alpha_group['pc'] /  alpha_group['pc_predicted']).mean()
    #     return pc_data['pc_predicted'] * c 
    pc_data = data.groupby(['alpha', 'n']).apply(get_pc).reset_index()
    pc_data.rename(columns={0: 'pc'}, inplace=True)
    pc_data['pc_predicted'] = pc_data['n'] ** (-2 * (1- pc_data['alpha']))
    pc_data['prediction_ratio'] = pc_data['pc'] /  pc_data['pc_predicted']
    pc_data['pc_predicted'] *= pc_data.groupby('alpha')['prediction_ratio'].transform('mean')
    # c = (pc_data['pc'] /  pc_data['pc_predicted']).mean()
    # pc_data['pc_predicted'] *= c
    return pc_data

