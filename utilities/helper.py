import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Union
from matplotlib.patches import FancyArrowPatch


BLUE, ORANGE, GREEN, RED, PURPLE, BROWN, PINK, GREY, KHAKI, AZURE = sns.color_palette()
YELLOW = (255/255, 192/255, 0/255)

CP = {	"blue" : BLUE,
	"orange" : ORANGE,
	"green" : GREEN,
	"red" : RED,
	"purple" : PURPLE,
	"brown" : BROWN,
	"pink" : PINK,
	"grey" : GREY,
	"khaki" : KHAKI,
	"azure" : AZURE,
    "yellow" : YELLOW}
 


def draw_network(
    graph: nx.Graph,
    sources: list=None,
    strc_color=BLUE,
    chain_color=GREY,
    source_color=RED,
    strc_size: float=300,
    chain_size: float=150,
    aspect_ratio: float=None,
    with_labels: Union[bool, str]="strc",
    arc_edges: bool=False,
    **args
):
    draw_args = {"G": graph}
    # define sources
    if sources is None:
        if "sources" in graph.graph:
            sources = graph.graph["sources"]
        else:
            sources = []
    # check for pos
    pos = nx.get_node_attributes(graph, "pos")
    if len(pos) != graph.number_of_nodes():
        pos = nx.kamada_kawai_layout(graph)
    draw_args["pos"] = pos
    # node_color
    if "node_color" not in args:
        node_color = [
            source_color if node in sources else strc_color if degree != 2 else chain_color
            for node, degree in graph.degree
        ]
        draw_args["node_color"] = node_color
    # node_size
    if "node_size" not in args:
        node_size = [
            strc_size if (node in sources) or (degree != 2) else chain_size
            for node, degree in graph.degree
        ]
        draw_args["node_size"] = node_size
    # labels
    if with_labels == "strc":
        labels={node: node for node, degree in graph.degree if degree != 2 or node in sources}
        draw_args["labels"] = labels
        draw_args["with_labels"] = True
        draw_args["font_color"] = 'white'
    elif isinstance(with_labels, bool):
        draw_args["with_labels"] = with_labels
    else:
        raise ValueError("with_labels should be bool or 'strc'")
    draw_args["alpha"] = 0.9
    # draw_args["linewidths"] = 1
    draw_args["edgecolors"] = "black"
    # ax
    if 'ax' not in args:
        args['ax'] = plt.gca()
    if aspect_ratio is not None:
        ax = args['ax'] if 'ax' in args else plt.gca()
        ax.set_aspect(aspect_ratio)
    draw_args.update(args)
    if arc_edges:
        cured_args = {
            'color': args.get("edge_color", "grey"),
            'rad': args.get("rad", 0.2),
            'width': args.get("width", 1.0),
            'alpha': args.get("alpha", 1.0),
        }
        draw_curved_edges(ax=draw_args['ax'], graph=graph, pos=pos, **cured_args)
        draw_args['width'] = 0
        if 'rad' in draw_args:
            del draw_args['rad']
        # draw_arc_edges(graph, pos, draw_args['ax'])
    nx.draw(**draw_args)


def draw_curved_edges(
    ax,
    graph,
    pos,
    *,
    rad=0.2,
    color="grey",
    width=1.0,
    alpha=0.8,
    zorder=0,
):
    """
    Draw curved edges using quadratic Bezier curves.
    """
    for u, v in graph.edges():
        (x1, y1) = pos[u]
        (x2, y2) = pos[v]

        patch = FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-",      # no arrow head
            linewidth=width,
            color=color,
            alpha=alpha,
            zorder=zorder,
        )
        ax.add_patch(patch)


