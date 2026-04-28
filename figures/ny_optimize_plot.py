import os
from typing import Optional, List, Tuple
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.axes import Axes

import optimal_network as ON
from figures.ny_plot import NEW_YORK_ASYMPTOTIC_FILE_NAME
import utilities.read_write as rw
from utilities.figures_utilities import (
    BLUE, cm_to_inch, configure_axes, format_p_scientific,
    net_plot, saidi_with_lengths,
)

# Constants
TEXT_SIZE = 17
LABEL_SIZE = 17
TICK_SIZE = 10
LINE_WIDTH = 2
MARKER_WIDTH = 1.5
LETTERS = ['a', 'b', 'c', 'd', 'e', 'f']


# ===== Plots =====

def ny_optimize_plot(nyc_optimal_path: str = r"data/nyc_optimal.nxjson", save: bool = False) -> plt.Figure:
    """
    Generate the full NY optimize plot grid.

    Args:
        nyc_optimal_path (str): Path to the optimal network simulation data.
        save (bool): Whether to save the figure to disk.

    Returns:
        plt.Figure: The generated figure.
    """
    fig, axs = plt.subplots(
        3, 2, figsize=(2 * 7 / cm_to_inch, 3 * 7 / cm_to_inch),
        gridspec_kw={'wspace': 0.4, 'hspace': 0.4}
    )
    
    # Read data
    data = rw.read_nxjson(nyc_optimal_path)
    data = data.explode(["p", "f"])  # Explode metrics
    pc_data = calculate_pc_data(data)

    # Plotting
    lineplots(pc_data, axs)
    draw_sample_graph(data, n_idx=7, axs=axs)
    number_plots(axs)

    # Save logic
    if save:
        base_save_dir = r'outputs\ny_optimize_plot'
        os.makedirs(base_save_dir, exist_ok=True)
        plt.savefig(
            os.path.join(base_save_dir, 'ny_optimize_plot.svg'),
            bbox_inches='tight',
            dpi=300,
            transparent=True
        )
        plt.savefig(
            os.path.join(base_save_dir, 'ny_optimize_plot.pdf'),
            bbox_inches='tight',
            dpi=300,
            transparent=True
        )
    return fig


def draw_sample_graph(data: pd.DataFrame, n_idx: int, axs: np.ndarray) -> None:
    """
    Draw network sample graphs on the provided axes grid.

    Args:
        data (pd.DataFrame): Simulation data.
        n_idx (int): Index for the desired graph size 'n'.
        axs (np.ndarray): 2D array of matplotlib Axes.
    """
    all_n = sorted(set(data['n']))
    all_alpha = sorted(set(data['alpha']))
    all_p = sorted(set(data['p']))
    
    n = all_n[n_idx]
    p = all_p[5]
    
    alpha_rows = {
        alpha: data.query('alpha == @alpha and n == @n and p == @p').iloc[0]
        for alpha in all_alpha
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


def lineplots(pc_data: pd.DataFrame, axs: np.ndarray) -> None:
    """
    Plot the critical probability transitions for different alphas.

    Args:
        pc_data (pd.DataFrame): Processed critical probability data.
        axs (np.ndarray): 2D array of matplotlib Axes.
    """
    for (alpha, group), ax in zip(pc_data.groupby('alpha'), axs[:, 1]):
        sns.lineplot(
            data=group, x='n', y='pc_predicted', color='black', linewidth=LINE_WIDTH, ax=ax
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


def number_plots(axs: np.ndarray) -> None:
    """
    Add reference letters (a, b, c, etc.) to the top left of each subplot.

    Args:
        axs (np.ndarray): Array of matplotlib Axes.
    """
    for i, ax in enumerate(axs.flatten()):
        ax.text(
            -0.03, 1.05, f"{LETTERS[i]}",
            fontsize=TEXT_SIZE, color="black",
            ha='right', va='bottom',
            transform=ax.transAxes,
        )


# ===== Simulations =====

def simulate_ny_optimal_networks(
    nyc_asymptotic_path: str = NEW_YORK_ASYMPTOTIC_FILE_NAME,
    nyc_optimal_path: str = r"data/nyc_optimal.nxjson",
    skip_existing_data: bool = True
) -> pd.DataFrame:
    """
    Simulate failures on optimal network designs for New York data.

    Args:
        nyc_asymptotic_path (str): Path to the base asymptotic NY data.
        nyc_optimal_path (str): Output path for the simulated optimal networks.
        skip_existing_data (bool): Whether to load existing data and resume.

    Returns:
        pd.DataFrame: DataFrame containing generated networks and their simulated SAIDI.
    """
    print("!!!TODO: set sources to optimal value!!!")
    alpha_list = [0.3, 0.4, 0.5]
    p_list = np.logspace(-4, -2, 10)
    df = rw.read_nxjson(nyc_asymptotic_path).query('trial == 0')
    
    if skip_existing_data and os.path.exists(nyc_optimal_path):
        data = rw.read_nxjson(nyc_optimal_path)
    else:
        data = pd.DataFrame(columns=["r", "alpha", "n", "p", "f", "graph", "sources"])
        
    data.set_index(['alpha', 'n'], inplace=True)
    
    for alpha in alpha_list:
        rng = np.random.default_rng(100)
        for _, row in df.iterrows():
            n = row.n
            graph = row.graph
            r = int(n ** alpha)
            
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
                sources=list(row.sources),
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
            print(f"Saved optimal NY simulation data to: {nyc_optimal_path}")
            rw.write_nxjson(data.reset_index(), nyc_optimal_path)
            
    return data


def calculate_pc_data(data: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate the critical probability transition data from the simulation outputs.

    Args:
        data (pd.DataFrame): Exploded simulation data containing 'p' and 'f'.

    Returns:
        pd.DataFrame: DataFrame containing empirical and predicted critical probabilities.
    """
    def get_pc(group: pd.DataFrame) -> float:
        group = group.sort_values('p')
        group = group.iloc[5:]
        c = np.mean(group['f'] / (group['p'] ** 2))
        return 1 / c

    pc_data = data.groupby(['alpha', 'n']).apply(get_pc).reset_index()
    pc_data.rename(columns={0: 'pc'}, inplace=True)
    
    # Calculate expected theoretical transition
    pc_data['pc_predicted'] = pc_data['n'] ** (-2 * (1 - pc_data['alpha']))
    
    # Align theoretical curve with empirical data vertically
    pc_data['prediction_ratio'] = pc_data['pc'] / pc_data['pc_predicted']
    pc_data['pc_predicted'] *= pc_data.groupby('alpha')['prediction_ratio'].transform('mean')
    return pc_data
