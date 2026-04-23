import seaborn as sns
import numpy as np
from typing import List, Tuple, Literal
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, AutoMinorLocator, LogLocator
from pyproj import Transformer


from utilities import draw_network
from indexes import GraphRel, edge_probs_by_length, RelSimulationResult
from indexes.probs import Float, Array
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import Voronoi
import contextily as ctx

def _norm(x: tuple) -> tuple:
    return tuple([a / 255 for a in x])

cm_to_inch = 2.54
RED =  sns.color_palette()[3] # _norm((165, 0, 33))
GREEN = sns.color_palette()[2] # _norm((145, 196, 110))
GREY = sns.color_palette()[7]
BLUE =  sns.color_palette()[0] # _norm((68, 114, 196))
YELLOW = _norm((255, 192, 0))

# Font size constants
TEXT_SIZE = 18
LABEL_SIZE = 18
TICK_SIZE = 14
ANNOTATION_SIZE = 16
LINE_WIDTH = 3

### axes ####

def configure_axes(ax, label_size=None, tick_size=None):
    if label_size is None:
        label_size = LABEL_SIZE    
    if tick_size is None:
        tick_size = TICK_SIZE

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(axis="both", which="major", labelsize=tick_size)
    ax.tick_params(axis="both", which="minor", length=4, color="grey", labelcolor="grey")

    ax.set_xlabel(ax.get_xlabel(), fontsize=label_size)
    ax.set_ylabel(ax.get_ylabel(), fontsize=label_size, rotation=90)

    # ---- tick logic ----
    for axis in (ax.xaxis, ax.yaxis):
        if axis.get_scale() == "linear":
            axis.set_minor_locator(AutoMinorLocator(2))  # exactly one minor tick
        else:
            # log / symlog / etc: let matplotlib or user decide
            axis.set_minor_locator(LogLocator(base=10, subs=range(2, 10)))


def format_p_scientific(p: float) -> str:
    # Zero
    if p == 0 or np.isclose(p, 0.0):
        return r"$0$"
    # Sign and magnitude
    sign = "-" if p < 0 else ""
    x = abs(p)
    # Scientific parts
    exponent = int(np.floor(np.log10(x)))
    coeff = x / (10 ** exponent)        # in [1, 10)
    # Round to 1 significant digit
    c = int(np.round(coeff))
    if c == 10:
        # rollover, e.g. 9.99 -> 10 -> 1 * 10^(e+1)
        c = 1
        exponent += 1

    # Build LaTeX: omit "1·" and keep sign
    coeff_str = f"{sign}{c}\\times " if c != 1 else sign
    return f"${coeff_str}10^{{{exponent}}}$"

def get_xylim_from_graphs(graphs: List[nx.Graph]) -> Tuple[int]:
    # Calculate bounding box for all graphs in EPSG:4326 coordinates
    all_x_coords = []
    all_y_coords = []
    for graph in graphs:
        # Get node positions from the graph
        pos = nx.get_node_attributes(graph, "pos")
        if pos:
            coords = np.array(list(pos.values()))
            all_x_coords.extend(coords[:, 0])
            all_y_coords.extend(coords[:, 1])
    
    # Calculate global xlim and ylim in EPSG:4326 coordinates
    if all_x_coords and all_y_coords:
        xlim = (min(all_x_coords), max(all_x_coords))
        ylim = (min(all_y_coords), max(all_y_coords))
    else:
        # Fallback if no positions found
        xlim = None
        ylim = None
    return xlim, ylim


def set_custom_scientific_format(ax, axis='x', factor=-5, offset_fontsize=TICK_SIZE, int_labels=True):
    def sci_fmt(x, pos):
        if int_labels:
            return f"{int(x/10**factor)}"
        return f"{x/10**factor:.2g}"
    if axis in ['x', 'both']:
        ax.xaxis.set_major_formatter(FuncFormatter(sci_fmt))
        ax.annotate(
            f"×10$^{{{factor}}}$",
            xy=(1, 0), xycoords=('axes fraction', 'axes fraction'),
            xytext=(0, -20), textcoords='offset points',
            ha='right', va='top',
            color='grey',
            fontsize=offset_fontsize
        )
    if axis in ['y', 'both']:
        ax.yaxis.set_major_formatter(FuncFormatter(sci_fmt))
        ax.annotate(
            f"×10$^{{{factor}}}$",
            xy=(0, 1), xycoords=('axes fraction', 'axes fraction'),
            xytext=(0, 20), textcoords='offset points',
            ha='left', va='top',
            color='grey',
            fontsize=offset_fontsize
        )


###### plot networks on map ####

def net_plot(
        graph: nx.Graph,
        ax: plt.Axes,
        sources: list | None=None,
        pad: float | list[float] = 0.05,
        to_3857: bool=True,
        basemap: bool = True
    ):
    """Plot a network on `ax` (projected to Web-Mercator) with optional basemap.

    Parameters
    ----------
    graph : nx.Graph
        Graph with node attribute `'pos'` in EPSG:4326 (lon, lat).
    ax : matplotlib.axes.Axes
        Axes to draw on.
    sources : iterable
        Nodes considered sources; these are highlighted.
    pad : float or sequence of 4 floats, optional
        Padding around the data used to set axis limits. If a single float `f`
        is provided it is interpreted as a symmetric fraction of the data
        range (applied equally to all four sides). If a 4-tuple/list is
        provided it must be in the order ``[north, east, south, west]``; each
        entry is a fraction of the corresponding axis range (e.g., north is a
        fraction of the y-range added to the top).
    basemap : bool
        Whether to add a background basemap (CartoDB Positron).
    """
    from utilities import draw_network


    # sources
    if (sources is None):
        if "sources" in graph.graph:
            sources = graph.graph["sources"]
        else:
            sources = []
    # --- Ensure positions are in 3857 for correct basemap + geometry ---
    pos = nx.get_node_attributes(graph, "pos")  # assumed (lon, lat)
    pos_3857 = _pos_4326_to_3857(pos) if to_3857 else pos

    # Draw network using the projected positions
    draw_network(
        graph,
        pos=pos_3857,
        node_color=[RED if node in sources else "grey" for node in graph.nodes],
        node_size=[20 if node in sources else 0 for node in graph.nodes],
        edge_color="grey",
        edgecolors=None,
        with_labels=False,
        ax=ax,
    )

    # Redundant edges
    T_edges = set(map(_key, _minimum_spanning_tree(graph).edges))
    re_edges = [e for e in graph.edges if _key(e) not in T_edges]
    draw_edges(
        graph,
        pos=pos_3857,
        edgelist=re_edges,
        width=1,
        edge_color=[BLUE] * len(re_edges),
        ax=ax,
        zorder=100,
    )

    # Limits (these must be in 3857 too if you set them!)
    _set_xy_limits_from_pos(ax, pos_3857, pad=pad)

    if basemap:
        ctx.add_basemap(ax, crs=3857, source=ctx.providers.CartoDB.Positron)

    # remove copyrights
    for text in ax.texts:
        if getattr(text, "get_text", None) and text.get_text().startswith("(C)"):
            text.set_text("")


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

def _pos_4326_to_3857(pos_dict):
    # pos_dict: {node: (lon, lat)}
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    out = {}
    for n, (x, y) in pos_dict.items():
        X, Y = transformer.transform(x, y)
        out[n] = (X, Y)
    return out

def _set_xy_limits_from_pos(ax, pos, pad=0.1, min_range=None):
    """Set x/y limits from node positions with flexible padding.

    Parameters
    ----------
    ax : matplotlib Axes
        Axes to set limits on.
    pos : dict
        Mapping node -> (x, y) coordinates.
    pad : float or sequence of 4 floats, default=0.1
        If a scalar, interpreted as a symmetric fraction of the data range
        applied equally to all sides (same behavior as before).
        If a sequence (north, east, south, west), each entry is a fraction of
        the corresponding axis range and applied directionally:
            - north: fraction of y-range added to ymax
            - east:  fraction of x-range added to xmax
            - south: fraction of y-range subtracted from ymin
            - west:  fraction of x-range subtracted from xmin
    min_range : float or None
        Minimal span (in data units) enforced on both axes.
    """
    if not pos:
        return

    coords = np.asarray(list(pos.values()))
    x = coords[:, 0]
    y = coords[:, 1]

    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())

    dx = xmax - xmin
    dy = ymax - ymin

    if min_range is not None:
        dx = max(dx, min_range)
        dy = max(dy, min_range)

    # Directional padding: accept a scalar or a 4-tuple/list/array
    if isinstance(pad, (list, tuple, np.ndarray)):
        if len(pad) != 4:
            raise ValueError("pad must be a scalar or a 4-length sequence [north, east, south, west]")
        north, east, south, west = pad
        # Convert fractional pads into absolute offsets
        padx_left = float(west) * dx
        padx_right = float(east) * dx
        pady_bottom = float(south) * dy
        pady_top = float(north) * dy
        ax.set_xlim(xmin - padx_left, xmax + padx_right)
        ax.set_ylim(ymin - pady_bottom, ymax + pady_top)
    else:
        # scalar symmetric padding
        padx = float(pad) * dx
        pady = float(pad) * dy
        print(f"{xmin, padx, xmax, padx=}")
        ax.set_xlim(xmin - padx, xmax + padx)
        ax.set_ylim(ymin - pady, ymax + pady)



def plot_graph_with_basemap(
        graph,
        crs='EPSG:4326', 
        ax=None,
        xlim=None,
        ylim=None,
        **draw_args
    ):
    import geopandas as gpd
    from pyproj import transformer
    if ax is None:
        ax = plt.gca()
    # Initialize a transformer to convert to Web Mercator
    pos = nx.get_node_attributes(graph, "pos").values()
    pos_geo = gpd.GeoSeries(gpd.points_from_xy(
        [p[0] for p in pos], [p[1] for p in pos]
    ), crs=crs).to_crs(4326)
    pos = {node: (p.x, p.y) for node, p in zip(graph.nodes, pos_geo)}
    # Draw the graph using the transformed positions
    draw_network(graph, ax=ax, pos=pos, **draw_args)
    # xy lims
    calculated_xlim = ax.get_xlim()
    calculated_ylim = ax.get_ylim()
    # Transform xlim and ylim if provided and if CRS is different from EPSG:4326
    if xlim is not None and crs != 'EPSG:4326':
        # Transform the xlim coordinates from EPSG:4326 to target CRS
        trans_points = gpd.GeoSeries(gpd.points_from_xy(xlim, [0, 0]), crs=crs).to_crs(4326)
        xlim = [p.x for p in trans_points]
        # xlim_min, xlim_max = xlim
        # xlim_min_transformed, _ = transformer.transform(xlim_min, 0)  # y doesn't matter for x transformation
        # xlim_max_transformed, _ = transformer.transform(xlim_max, 0)
        # xlim = (xlim_min_transformed, xlim_max_transformed)
    elif xlim is None:
        xlim = calculated_xlim
        
    if ylim is not None and crs != 'EPSG:4326':
        # Transform the ylim coordinates from EPSG:4326 to target CRS
        trans_points = gpd.GeoSeries(gpd.points_from_xy([0, 0], ylim), crs=crs).to_crs(4326)
        ylim = [p.y for p in trans_points]
        # ylim_min, ylim_max = ylim
        # _, ylim_min_transformed = transformer.transform(0, ylim_min)  # x doesn't matter for y transformation
        # _, ylim_max_transformed = transformer.transform(0, ylim_max)
        # ylim = (ylim_min_transformed, ylim_max_transformed)
    elif ylim is None:
        ylim = calculated_ylim
        
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    # Add the basemap
    # ctx.add_basemap(ax, crs=crs, source=ctx.providers.OpenStreetMap.Mapnik)
    ctx.add_basemap(ax, crs=crs, source=ctx.providers.CartoDB.Positron)
    # ctx.add_basemap(ax, crs=crs, source=ctx.providers.CartoDB.VoyagerNoLabels)
    # ctx.add_basemap(ax, crs=crs, source=ctx.providers.Esri.WorldTerrain)

    for text in ax.texts:
        text._text = ""


#### SAIDI ###
def saidi_with_lengths(
        G: nx.Graph,
        sources: list,
        p: float | list[float],
        mode: str = "rate",
        T_days: float = 365 * 5,
        mean_cycle_days: float = 0.5,
        rng: np.random.Generator = None,
        show_progress: bool = True,
    ) -> float | list[float]:
    """Compute SAIDI (time-average disconnected fraction) for a graph where edge
    failure probabilities are derived from geometric length.

    Behavior
    - If `p` is a single scalar: treat it as a rate (if mode='rate') or a target
      mean failure probability (if mode='mean') and return one float (mean SAIDI).
    - If `p` is an iterable: compute SAIDI for each value and return a list of
      floats in the same order.

    Implementation notes
    - Uses `edge_probs_by_length(..., pos_attr='pos', prob_attr=None, tol=1e-6)` to
      compute per-edge failure probabilities.
    - Wraps probabilities with `Float` (scalar case) or `Array` (vector case) and
      uses `GraphRel.calc_rel_simulation` to run the SAIDI simulations.

    Parameters
    - G: networkx.Graph with node attribute `'pos'` = (x, y) for every node.
    - sources: list of source nodes (passed directly to GraphRel).
    - p: float or iterable of floats (see mode description).
    - mode: 'rate' or 'mean' (see `edge_probs_by_length`).
    - T_days, mean_cycle_days: simulation timing parameters forwarded to simulator.
    - rng: optional numpy Generator for reproducibility.
    - show_progress: show progress bars while simulating when True.

    Returns
    - float if a single p was supplied, otherwise list[float].
    """
    # --- sanity checks ---
    if mode not in ("rate", "mean"):
        raise ValueError("mode must be 'rate' or 'mean'")

    # Normalize p to a list for uniform processing
    is_scalar = np.isscalar(p)
    p_list = [float(p)] if is_scalar else [float(x) for x in p]
    if len(p_list) == 0:
        raise ValueError("p must contain at least one value")

    # Compute edge probabilities for each p value using the requested settings.
    # edge_probs_by_length returns (probs_dict, rate_used)
    probs_and_rates = [edge_probs_by_length(G, p=val, prob_attr=None, pos_attr="pos", mode=mode, tol=1e-7) for val in p_list]
    probs_list = [pr[0] for pr in probs_and_rates]

    # Ensure consistent edge keys across all computations
    base_keys = set(probs_list[0].keys())
    for i, pr in enumerate(probs_list[1:], start=1):
        if set(pr.keys()) != base_keys:
            raise ValueError(f"Edge set mismatch between p values at index {i}")

    # Build GraphRel and run simulation
    if is_scalar:
        # Scalar case: wrap each prob in Float and run one simulation
        edges_prob = {e: Float(float(pe)) for e, pe in probs_list[0].items()}
        gr = GraphRel(G, sources=sources, edges_prob=edges_prob)
        res = gr.calc_rel_simulation(T_days=T_days, mean_cycle_days=mean_cycle_days, rng=rng, show_progress=show_progress)
        return float(res.rel_result)

    # Vector case: build Array objects per-edge where each element corresponds to a p
    edges_prob = {e: Array([pr[e] for pr in probs_list]) for e in probs_list[0].keys()}
    gr = GraphRel(G, sources=sources, edges_prob=edges_prob)
    res = gr.calc_rel_simulation(T_days=T_days, mean_cycle_days=mean_cycle_days, rng=rng, show_progress=show_progress)

    # `calc_rel_simulation` returns a tuple of SaidiSimulationResult in the Array case
    return [float(r.rel_result) for r in res]

def failing_rate_from_spanning_tree(points: np.array, p_mean: float) -> float:
    tree = _euclidean_spanning_tree(points)
    _, rate = edge_probs_by_length(tree, p=p_mean, pos_attr="pos", mode="mean", tol=1e-6)
    return rate


def _euclidean_spanning_tree(points: np.ndarray) -> nx.Graph:
    """
    Return the euclidean minimal spanning tree of a sef of points
    use the fact that the EMST is a subgaph of the delaunay_graph.

    The Return graph has a node attribute "pos" and weight attribute weight
    """
    points = np.asarray(points, dtype="float")
    # build Delaunay graph
    vor = Voronoi(points)
    edges = np.array(vor.ridge_points, dtype="int")
    dists = np.linalg.norm(
        points[edges[:, 0], :] - points[edges[:, 1], :],
        axis=1
    )
    delaunay_graph = nx.Graph()
    delaunay_graph.add_weighted_edges_from([(*e, d) for e, d in zip(edges, dists)])
    # find spanning tree
    T = nx.minimum_spanning_tree(delaunay_graph, weight="weight")
    T = nx.edge_subgraph(delaunay_graph, T.edges)
    nx.set_node_attributes(T, dict(enumerate(points)), "pos")
    return T

from pathlib import Path

def save_axs_without_text(axs: List[plt.Axes], path: str | Path):
    """
    Takes a list of axes, removes all text except ticks, 
    and saves each axis as a separate SVG in the given path.
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)

    for i, ax in enumerate(axs):
        # Remove titles and labels
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel("")
        
        # Remove any other text or annotations
        for txt in ax.texts:
            txt.set_visible(False)
            
        # Remove legend if it exists
        legend = ax.get_legend()
        if legend is not None:
            legend.set_visible(False)

        # Get the bounding box of the axis including ticks and tick labels
        fig = ax.figure
        
        # Hide all other axes to prevent them from being embedded in the SVG
        visibility_map = {}
        for other_ax in fig.axes:
            visibility_map[other_ax] = other_ax.get_visible()
            if other_ax != ax:
                other_ax.set_visible(False)

        # In some backend situations, we need to draw the figure first to get the renderer
        if getattr(fig.canvas, "get_renderer", None) is not None:
            renderer = fig.canvas.get_renderer()
        else:
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            
        extent = ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted())
        
        # Save as SVG
        file_path = path / f"ax_{i}.svg"
        fig.savefig(file_path, bbox_inches=extent, format='svg', transparent=True)
        
        # Restore visibility
        for other_ax in fig.axes:
            other_ax.set_visible(visibility_map[other_ax])

                    