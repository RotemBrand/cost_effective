import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.spatial import ConvexHull, voronoi_plot_2d, Voronoi
from matplotlib.patches import Polygon

from indexes.utilities import create_sparse_graph_with_pos
from utilities.helper import draw_network
import optimal_network as ON
from utilities.figures_utilities import BLUE, GREEN, GREY, RED, cm_to_inch, TEXT_SIZE
import optimal_network.construct_strc as CS
from indexes.utilities import get_skeleton_graph
from figures.optimize_plot import plot_line_with_title

NODES_LINEWIDTH = 1.5
STRC_SIZE = 100
WIDTH = 1.5


def example_draw(save: bool=False):
    # figure
    fig, axs = plt.subplots()
    graph = build_example_graph()

    # draw
    draw_network(graph, chain_size=0, strc_size=90, with_labels=False)
    if save:
        file_name = r'outputs\supp_graph_example.svg'
        plt.savefig(
            file_name,
            bbox_inches='tight',
            dpi=300
        )
        print(f"graph saved to {r'outputs\supp_graph_example.svg'}")


def build_example_graph():
    """Build an example network to show how to calculate reliability"""
    # build strc
    strc = nx.from_edgelist(
        list(nx.complete_graph(4).edges) +
        [(1, 4)] +
        nodes_from([(0, 1), (1, 2), (0, 3), (3, 2), (0, 4), (4, 2)], start=4) + 
        nodes_from(nx.random_tree(n=7, seed=10).edges, start=8, mapping={0: 3}) +
        # [(2, 6)] + 
        nodes_from(nx.complete_graph(4).edges, start=15) + 
        [(22, 6)] + 
        nodes_from(nx.complete_graph(4).edges, start=19) +
        [(19, 0), (22, 16)]
    )
    strc.remove_edges_from([(3, 0)])
    # strc = nx.from_edgelist([
    #     (0, 1), (0, 2), (1, 2), (2, 6),
    #     (1, 3), (2, 3), (3, 5), (5, 10), (10, 11), (11, 5),
    #     # (3, 4)
    # ])
    # add chains
    pos = nx.kamada_kawai_layout(strc)
    nx.set_node_attributes(strc, pos, "pos")
    rng = np.random.default_rng(10)
    strc.nodes[5]["pos"] = (strc.nodes[8]["pos"] + strc.nodes[7]["pos"]) / 2
    graph = create_sparse_graph_with_pos(
        strc,
        {e: int(5 * rng.random()) for e in strc.edges}
    )
    return graph



def nodes_from(edges, start=0, mapping=None):
    if not mapping:
        mapping = {}
    def change(v):
        if v in mapping:
            return mapping[v]
        return v + start
    return [[change(v) for v in e] for e in edges]


def build_example_graph2():
    strc = nx.from_edgelist([
        (0, 1), (0, 2), (2, 6), (1, 3), (3, 4), (3, 5), (5, 5)   
    ])



########## algorith, example ######

def supp_plot_spatial_network(save: bool=False):
    fig, axs = plt.subplots(
        3, 3, figsize=(25 / cm_to_inch, 30 / cm_to_inch),
        gridspec_kw={'hspace': 0.4}
    )
    # clusters
    plot_line_with_title(axs[0, 0], "Chains from points", color='black')
    axs[0, 0].set_title("Points")
    _draw_building_method(draw_points=True, ax=axs[0, 0], points_color=GREY)
    axs[0, 1].set_title("Random 3R centers")
    _draw_building_method(draw_points=True, draw_cluster_centers=True, kmeans_max_iter=1, ax=axs[0, 1], points_color=GREY)
    axs[0, 2].set_title(r"Assign equal clusters")
    _draw_building_method(draw_points=True, draw_cluster_centers=True, draw_clusters=True, kmeans_max_iter=7, ax=axs[0, 2])


    # strc
    plot_line_with_title(axs[1, 0], "Structure graph", color='black')
    axs[1, 0].set_title("Voronoi")
    _draw_building_method(draw_points=False, draw_cluster_centers=True, draw_voronoi=True, kmeans_max_iter=7, ax=axs[1, 0])
    axs[1, 1].set_title("Trips")
    _draw_building_method(draw_points=False, draw_cluster_centers=True, draw_trips=True, kmeans_max_iter=7, ax=axs[1, 1])
    axs[1, 2].set_title("Structure graph")
    _draw_building_method(draw_points=False, draw_cluster_centers=True, draw_strc=True, kmeans_max_iter=7, ax=axs[1, 2])

    # chains
    plot_line_with_title(axs[2, 0], "Connect points", color='black')
    axs[2, 0].set_title("Points to chains")
    _draw_building_method(draw_points=True, draw_clusters=True, draw_some_chains=True, kmeans_max_iter=7, ax=axs[2, 0])
    # _draw_building_method(draw_points=True, draw_strc=False, draw_clusters=True, draw_graph=True, kmeans_max_iter=7, ax=axs[2, 0])
    axs[2, 1].set_title("Optimal network")
    _draw_building_method(draw_points=True, draw_graph=True, kmeans_max_iter=7, ax=axs[2, 1])
    axs[2, 2].set_title("Local improvements")
    _draw_building_method(draw_points=True, draw_improvment=True, kmeans_max_iter=7, ax=axs[2, 2])

    # configure axes
    for ax in axs.flatten():
        ax.axis(False)
    number_plots(axs)
    # save
    if save:
        fig.savefig(
            r'outputs\supp_spatial_network.svg',
            bbox_inches='tight',
            dpi=300,
            transparent=True
        )



def _draw_building_method(
        draw_points: bool=True,
        draw_clusters: bool=False,
        draw_cluster_centers: bool=False,
        draw_voronoi: bool=False,
        draw_trips: bool=False,
        draw_strc: bool=False,
        draw_graph: bool=False,
        draw_some_chains: bool=False,
        draw_improvment: bool=False,
        kmeans_max_iter: int=7,
        points_color=GREEN,
        ax=None,
    ):
    if ax is None:
        ax = plt.gca()
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
        points, m, max_iter=kmeans_max_iter, chain_len_sigma=0, random_state=seed
    )

    data = pd.DataFrame(np.c_[points, chains], columns=["x", "y", "c"])
    data['node'] = list(range(len(points)))
    data['node'] = data['node'].apply(int).astype('int')
    data['c'] = data['c'].astype(int)

    # pick an ordered palette and bind each cluster id -> color
    clusters = list(set(data['c']))
    color_map = {c: GREEN for c in clusters}
    data['color'] = data['c'].map(color_map)

    # scatter with the exact same palette mapping (legend hidden as before)
    if draw_points:
        sns.scatterplot(
            x=data['x'], y=data['y'], ax=ax,
            c='white', edgecolors=points_color,
            s=30, zorder=1, linewidths=NODES_LINEWIDTH
        )
    
    if draw_cluster_centers:
        sns.scatterplot(
            x=centers[:, 0], y=centers[:, 1], ax=ax,
            s=300, c='black', marker='X',
        )


    if draw_clusters:
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
        trips_G.add_edge(s, str(c))
        trips_G.add_edge(str(c), t)
        trips_G.nodes[s]["pos"] = points[s]
        trips_G.nodes[t]["pos"] = points[t]
        trips_G.nodes[str(c)]["pos"] = centers[c]

    if draw_voronoi:
        vor = Voronoi(centers)
        voronoi_plot_2d(vor, ax=ax, show_points=False)
        vertex_to_point = _get_closest_point_to_each_vertex(points=points, centers=centers)
        selection_G = nx.from_edgelist(vertex_to_point, create_using=nx.DiGraph)
        nx.draw_networkx_edges(
            selection_G,
            pos={v: v for v in selection_G.nodes},
            node_size=0,
            edge_color="black",
            ax=ax
        )
        selected_points = np.array([np.array(point) for _, point in vertex_to_point])
        sns.scatterplot(
            x=selected_points[:, 0],
            y=selected_points[:, 1],
            s=50,
            c='white', edgecolor=GREY, linewidths=NODES_LINEWIDTH,
            ax=ax
        )
        sns.scatterplot(
            x=vor.vertices[:, 0],
            y=vor.vertices[:, 1],
            c='red',
            s=50,
            zorder=100,
            edgecolor='black',
            ax=ax
        )

    if draw_trips:
        trips_graph = nx.from_edgelist(chosen_trips, create_using=nx.MultiGraph)
        matching = nx.maximal_matching(nx.Graph(trips_graph))
        for s, t in matching:
            c = next(iter(trips_graph[s][t]))
            data = np.array([points[s], centers[c], points[t]])
            sns.lineplot(
                x=data[:2, 0] , y=data[:2, 1], c=GREY, linewidth=2, ax=ax
            )

            sns.lineplot(
                x=data[1:, 0] , y=data[1:, 1], c=GREY, linewidth=2, ax=ax
            )
            sns.scatterplot(
                x=data[[0, -1], 0], y=data[[0, -1], 1], s=100,
                color='white', edgecolor=BLUE, zorder=100, ax=ax, linewidths=NODES_LINEWIDTH    
            )


    if draw_strc:
        draw_network(
            trips_G,
            node_color=['white' if d != 2 else GREY for node, d in trips_G.degree],
            edgecolors=[GREY if d == 2 else BLUE for node, d in trips_G.degree],
            edge_color='grey',
            linewidths=NODES_LINEWIDTH,
            chain_size=0,
            strc_size=STRC_SIZE,
            with_labels=False,
            width=WIDTH + 1,
            ax=ax,
        )
    
    ######## graph ########
    if draw_some_chains:
        subgraph = ON.add_chains_to_strc(points, chains, chosen_trips, max_init_iter=2, debug=False)
        strc = get_skeleton_graph(subgraph, sources=[0])
        _set_weights_to_strc(strc)
        # choose chains
        chains = list(nx.min_weight_matching(nx.Graph(strc), weight="weight"))[:2]
        edges = []
        for chain in chains:
            edges += list(map(
                lambda x: x[:2],
                strc.edges[(*chain, 0)]["subgraph"].edges
            ))

        subgraph = nx.edge_subgraph(subgraph, edges)
        draw_network(
            subgraph,
            node_color='white',
            edgecolors=[GREEN if d == 2 else BLUE for node, d in subgraph.degree],
            strc_size=STRC_SIZE,
            chain_size=30,
            with_labels=False,
            edge_color='grey',
            linewidths=NODES_LINEWIDTH,
            width=WIDTH,
            ax=ax
        )

    if draw_graph or draw_improvment:
        subgraph = ON.add_chains_to_strc(points, chains, chosen_trips, max_init_iter=2, debug=False)
    if draw_graph:
        draw_network(
            subgraph,
            node_color='white',
            strc_size=STRC_SIZE,
            chain_size=30,
            with_labels=False,
            edge_color='grey',
            linewidths=NODES_LINEWIDTH,
            edgecolors=[GREEN if d == 2 else BLUE for node, d in subgraph.degree],
            width=WIDTH,
            ax=ax
        )
        return subgraph
    if draw_improvment:
        from optimal_network import optimize_rel_weight_ratio
        new_subgraph = optimize_rel_weight_ratio(subgraph, max_risk_gain=0.1, max_changes=5)
        edges_diff = _get_edges_diff(subgraph, new_subgraph)
        draw_network(
            new_subgraph,
            node_color='white',
            strc_size=STRC_SIZE,
            chain_size=30,
            with_labels=False,
            edge_color=[RED if tuple(sorted(e)) in edges_diff else 'grey' for e in new_subgraph.edges],
            linewidths=NODES_LINEWIDTH,
            edgecolors=[GREEN if d == 2 else BLUE for node, d in new_subgraph.degree],
            width=WIDTH,
            ax=ax
        )

    return None


def _set_weights_to_strc(strc: nx.MultiGraph):
    """Set weight attribute for each chain in the strc by the total len of its subgraph"""
    for chain, data in strc.edges.items():
        subgraph = data["subgraph"]
        weight = 0
        for e in subgraph.edges:
            weight += np.linalg.norm(subgraph.nodes[e[0]]["pos"] - subgraph.nodes[e[1]]["pos"])
        strc.edges[chain]["weight"] = weight

def _get_edges_diff(graph1: nx.Graph, graph2: nx.Graph) -> set:
    """Return all the edges that are in graph2 and not in graph1"""
    edges1 = set([tuple(sorted(e)) for e in graph1.edges])
    edges2 = set([tuple(sorted(e)) for e in graph2.edges])
    return edges2 - edges1

def _get_closest_point_to_each_vertex(points: np.array, centers: np.array):
    vor = Voronoi(centers)
    centers_dict = dict(enumerate(centers))

    # filter only inner vertices
    vertices = vor.vertices
    vertices = vertices[vertices[:, 0] <= centers[:, 0].max()]
    vertices = vertices[vertices[:, 0] >= centers[:, 0].min()]
    vertices = vertices[vertices[:, 1] <= centers[:, 1].max()]
    vertices = vertices[vertices[:, 1] >= centers[:, 1].min()]
    chosen_points_id = CS._closest_point_to_each_vertex(points=points, vertices=vertices)
    vertex_to_point = tuple(zip(
        map(tuple, vertices),
        map(tuple, points[chosen_points_id])
    ))
    return vertex_to_point

LETTERS = ['a', 'b', 'c' ,'d', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']
def number_plots(axs):
    for i, ax in enumerate(axs.flatten()):
        ax.text(
            -0.03, 1.05, f"{LETTERS[i]}",
            fontsize=TEXT_SIZE, color="black",
            ha='right', va='bottom',
            transform=ax.transAxes,
        )

