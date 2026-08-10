import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from epyt import epanet
from typing import Dict
import os
import pandas as pd
import pandapower as pp
import geopandas as gpd
from shapely.geometry import LineString
import json
from shapely.geometry import shape
from indexes import GraphRel, get_skeleton_graph, defult_sources, edge_probs_by_length
from utilities import draw_network
from indexes.simulation import RelType
from optimal_network import optimal_network_from_points, optimize_rel_weight_ratio, improve_tree
from scipy.spatial.distance import pdist, squareform
import topohub
import utilities.read_write as rw

DATA_PATH = r"data/real_networks.nxjson"

#### plots #####

    
    


##### main #####

def main_get_real_networks_and_improve(file_name: str=DATA_PATH) -> pd.DataFrame:
    """Read real-world networks, compute improved optimal networks and save results.

    This is a convenience top-level function that gathers networks from various
    sources (TopoHub, water, power), runs optimization routines to compute an
    optimal network for each, and optionally saves the resulting DataFrame to
    `file_name` using the project's read/write utilities.

    Parameters
    ----------
    file_name : str
        Path to write the resulting DataFrame as JSON (passed to RW.pd_to_json).

    Returns
    -------
    pd.DataFrame
        DataFrame with original and optimized networks and derived metrics.
    """
    print("==== read networks ======")
    networks_df = get_all_relevent_networks()
    print("==== optimize networks ======")
    networks_df = get_optimal_networks(networks_df)


    # save
    networks_df = _remove_unjson_types(networks_df)
    if file_name:
        rw.write_nxjson(networks_df, file_name)
    return networks_df

def _remove_unjson_types(networks_df: pd.DataFrame) -> pd.DataFrame:
    def fix_graph(G: nx.Graph) -> nx.Graph:
        G_copy = G.copy()
        for node in G_copy.nodes:
            G_copy.nodes[node]["pos"] = tuple(G_copy.nodes[node]["pos"])
        for e in G_copy.edges:
            if 'geometry' in G_copy.edges[e]:
                del G_copy.edges[e]["geometry"]
        return G_copy
    networks_df[["graph", "optimal_network"]] = networks_df[["graph", "optimal_network"]].map(fix_graph)
    return networks_df

def get_all_relevent_networks():
    """Collects networks from multiple sources into a single DataFrame.

    Searches topological datasets, water distribution examples and power
    distribution examples, converts them to NetworkX graphs and returns a
    concatenated DataFrame with metadata for each network.
    """
    networks_data = []
    # topohub
    topo_dict = {
        'Wavenet': nx.node_link_graph(topohub.get(f"topozoo/VtlWavenet2011")),
        'NetworkUSA': nx.node_link_graph(topohub.get(f"topozoo/NetworkUsa")),
    }
    topo_df = networks_dict_to_df(topo_dict)
    topo_df['type'] = 'communication'

    # water
    water_dict = get_all_water_networks(["Balerma", "Rural"])
    water_df = networks_dict_to_df(water_dict)
    water_df['type'] = 'water'

    # power
    power_dict = get_sfo_networks(improve=True)
    power_df = networks_dict_to_df(power_dict)
    power_df['type'] = 'power'

    # merge
    data = pd.concat([topo_df, water_df, power_df])
    return data

def get_optimal_networks(networks_df :pd.DataFrame) -> pd.DataFrame:
    """Compute and attach optimized networks and summary metrics to a DataFrame.

    Parameters
    ----------
    networks_df : pd.DataFrame
        DataFrame with at least columns ['name', 'type', 'graph'] describing
        the source networks. This function will add several columns describing
        the optimized networks and derived metrics.

    Returns
    -------
    pd.DataFrame
        Input DataFrame augmented with optimization results.
    """

    def get_optimal_network_by_type(row: pd.Series):
        print(f"===== {row['name']} =====")
        if row["type"] == "communication":
            oG = rw.dict_read_json(rf"data/zoo/optimal_{row['name']}.nxjson")
        elif row["type"] == "water":
            oG = get_optimal_network(row.graph, R=5)
            if row["name"] == "Balerma":
                oG = optimize_rel_weight_ratio(oG, max_risk_gain=0.5, max_changes=300)
        elif row["type"] == "power":
            if row["name"] == "SFO Davidson":
                oG = get_optimal_network(row.graph, R=6, debug=True, strc_n_init_iters=6, chain_n_init_iters=6, seed=10)
                oG = optimize_rel_weight_ratio(oG, max_risk_gain=0.3, max_changes=300)
            else:
                oG = get_optimal_network(row.graph, R=8, strc_n_init_iters=5, chain_n_init_iters=6, debug=True)
                oG = optimize_rel_weight_ratio(oG, max_risk_gain=1, max_changes=300)
        else:
            raise ValueError("Invalid type")
        if "sources" in row.graph.graph:
            sources = [x for x in row.graph.graph["sources"] if x in oG.nodes]
            if len(sources) > 0:
                oG.graph["sources"] = sources
        return oG
    
    rng = np.random.default_rng(100)
    # optimal network    
    networks_df['optimal_network'] = networks_df.apply(get_optimal_network_by_type, axis=1)
    networks_df["R_optimal_network"] = networks_df["optimal_network"].apply(lambda G: len(G.edges) - len(G) + 1)
    networks_df["optimal_network_weight"] = networks_df['optimal_network'].apply(get_total_weight)

    # saidis
    networks_df["rel_type"] = networks_df["type"].map(lambda t: "pairwise" if t == "communication" else "saidi")
    saidis = networks_df.apply(
        lambda row: get_G_and_OG_saidi(row["graph"], row["optimal_network"], rng=rng, rel_type=row["rel_type"]),
        axis=1
    )
    networks_df["saidi"] = list(map(lambda x: x[0], saidis))
    networks_df["optimal_network_saidi"] = list(map(lambda x: x[1], saidis))

    # ratios
    networks_df["R_ratio"] = 1 - networks_df["R_optimal_network"] / networks_df["R"]
    networks_df["weight_ratio"] = 1 - networks_df["optimal_network_weight"] / networks_df["total_weight"]
    networks_df["rel_ratio"] = 1 - networks_df["optimal_network_saidi"] / networks_df["saidi"]
    return networks_df
        


############################
#          analyze         #
############################
def networks_dict_to_df(net_dict: Dict[str, nx.Graph]) -> pd.DataFrame:
    """Convert a dict of name->NetworkX graph into a summary DataFrame.

    The returned DataFrame contains precomputed metrics (N, M, R, total weight,
    SAIDI, etc.) and the original graph object in the 'graph' column.
    """
    # read data
    networks_data = []
    for name, graph in net_dict.items():
        # size
        N = len(graph)
        M = len(graph.edges)
        R = M - N + 1
        # chains
        optimal_c = (N - 2 * (R - 1)) / (3 * (R - 1)) if R > 1 else N
        weight = get_total_weight(graph)
        chains = get_ext_chains(nx.Graph(graph))
        prec2_conn = get_prec_2_conn(nx.Graph(graph))
        p = 5e-4
        n_sources = len(graph.graph["sources"]) if "sources" in graph.graph else 1 
        optimal_saidi = get_optimal_saidi(N, R, p, n_sources)
        networks_data.append([name, N, M, R, R/N, optimal_c, weight, optimal_saidi, prec2_conn, chains, graph])
    # make df
    networks_pd = pd.DataFrame(networks_data, columns=["name", "N", "M", "R", "R/N", "optimal_c", "total_weight", "optimal_saidi", "prec2_conn", "ext_chains_len", "graph"])
    networks_pd.sort_values("optimal_c", ascending=False, inplace=True)
    return networks_pd


DEFUALT_ATTR = ["N", "R", "R/N", "c", "prec2_conn"]
def explore_plot(networks_pd: pd.DataFrame):
    """Small convenience plotting helper to visualize a list of networks.

    Plots networks in a grid using :func:`draw_network` and annotates basic
    statistics in the subplot title.
    """
    n = len(networks_pd)
    n_rows = int(np.ceil(n / 3))
    fig, axs = plt.subplots(n_rows, 3, figsize=(5*3, 5*n_rows))

    for (_, row), ax in zip(networks_pd.iterrows(), axs.flatten()):
        draw_network(row.graph, with_labels=False, chain_size=5, strc_size=10, ax=ax, edgecolors=None)
        name, N, R, c, chains, prec2_conn, saidi, optimal_saidi = row[["name", "N", "R", "optimal_c", "ext_chains_len", "prec2_conn", "saidi", "optimal_saidi"]]
        saidi_ratio = saidi / optimal_saidi
        ax.set_title(f"{name}\n{N=}, {R=}, {R/N=:.0%}, {c=:.0f}, 2_conn={prec2_conn:.1%}\nsaidi_ratio={saidi_ratio:.1f}, chains: {chains[:5]}")


def get_optimal_network(graph: nx.Graph, R=None, c_mean=None, **args) -> nx.Graph:
    """Create an optimal network (approximation) for a given graph.

    This wraps :func:`optimal_network_from_points` converting the input graph's
    node positions to a point cloud and selecting sensible defaults for the
    optimization parameters.
    """
    N, L = len(graph), len(graph.edges)
    if (R is None) and (c_mean is None):
        R = L - N + 1
    elif (R is None) and (c_mean is not None):
        R = int(N / (3*c_mean + 2) + 1)
        R = max([R, 1])
    print(f"N={len(graph)}, {R=}")
    points = np.array(list(nx.get_node_attributes(graph, "pos").values()))
    final_args = dict(strc_n_init_iters=10, seed=5)
    final_args.update(args)
    optimal_network = optimal_network_from_points(points, r=R, **final_args)
    return optimal_network


def get_total_weight(G: nx.Graph) -> float:
    """Return total Euclidean length of graph edges using node 'pos' attributes.

    If positions are missing this will raise a KeyError.
    """
    pos = nx.get_node_attributes(G, "pos")
    weights = [
        np.linalg.norm(np.array(pos[e[0]]) - np.array(pos[e[1]]))
        for e in G.edges
    ]
    return sum(weights)




def get_ext_chains(graph: nx.Graph) -> list:
    """Return a sorted list of external chain lengths for 2-edge-connected components.

    Contracts components and extracts chain lengths from skeleton graphs.
    """
    comps2 = list(nx.k_edge_components(graph, 2))
    graph_copy = graph.copy()
    chains = []
    for comp in comps2:
        strc_comp = get_skeleton_graph(nx.subgraph(graph, comp))
        comp_chains = list(nx.get_edge_attributes(strc_comp, "length").values())
        chains.extend(comp_chains)
        comp = list(comp)
        for node in comp[1:]:
            nx.contracted_nodes(graph_copy, comp[0], node, self_loops=False, copy=False)
    strc = get_skeleton_graph(graph_copy)
    chains += list(nx.get_edge_attributes(strc, "length").values())
    return sorted(chains)[::-1]

def get_prec_2_conn(graph: nx.Graph) -> float:
    """Return the fraction of nodes that belong to a non-trivial 2-edge-connected component.

    The result is in [0,1].
    """
    comps2 = filter(lambda comp: len(comp) > 1, nx.k_edge_components(graph, 2))
    return sum(len(comp) for comp in comps2) / len(graph)

def get_chain_rel(n: int, p: float) -> float:
    """Analytic reliability formula for a chain of length n at failure probability p."""
    return 1 / 6 * p ** 2 * (n **3 + 3 * n**2 + 2*n)

def get_optimal_saidi(N: int, R: int, p: float, n_sources: int=1) -> float:
    """Estimate SAIDI for the idealized optimal network given parameters.

    Parameters
    ----------
    N : int
        Number of nodes in the system.
    R : int
        Redundancy parameter (number of cycles + 1).
    p : float
        Failure probability used in the analytic chain formula.
    n_sources : int
        Number of sources in the network.
    """
    n = N - n_sources
    if R > 1:
        n_strc = 2 * (R - 1)
        m_strc = 3 * (R - 1)  + n_sources - 1
    elif R == 1:
        n_strc = 1
        m_strc = 1
    else:
        n_strc = n
        m_strc = n - 1
    l, n_long_chains = divmod(n - n_strc, m_strc)
    optimal_rel = (
        get_chain_rel(l, p) * (m_strc - n_long_chains) +
        get_chain_rel(l+1, p) * n_long_chains
    ) / N
    return optimal_rel


def get_G_and_OG_saidi(G: nx.Graph, OG: nx.Graph, rng: np.random.Generator, rel_type: RelType) -> list[float, float]:
    from indexes.probs import Float
    # get all probs
    G_prob, p_rate = edge_probs_by_length(G, p=5e-4, mode="mean")  
    G_prob = {e: Float(p) for e, p in G_prob.items()}
    OG_prob, _ = edge_probs_by_length(OG, p=p_rate, mode="rate")
    OG_prob = {e: Float(p) for e, p in OG_prob.items()}

    # calculate saidi per graph
    graphs = [G, OG]
    edge_probs = [G_prob, OG_prob]
    rels = []
    for graph, edge_prob in zip(graphs, edge_probs):
        sources = graph.graph.get("sources", [defult_sources(graph)])
        gr = GraphRel(graph, edges_prob=edge_prob, sources=sources)
        res = gr.calc_rel_simulation(rel_type=rel_type, T_days=365*5, mean_cycle_days=0.01, rng=rng)
        rels.append(res.rel_result)
    return rels




############################
#          water           #
############################
def water_inp_to_nx(inp_path: str) -> nx.Graph:
    """
    Convert an EPANET .inp water network (via epyt) into a NetworkX graph.

    Each node includes:
        - type: 'junction', 'reservoir', or 'tank'
        - demand, elevation, coordinates
        - is_source (True for reservoirs/tanks)
        - is_consumer (True if demand > 0)
        - weight (by default = demand)
    
    Each edge includes:
        - type: 'pipe', 'pump', or 'valve'
        - from, to
        - length, diameter, roughness (if available)
        - weight (by default = length)

    Parameters
    ----------
    inp_path : str
        Path to the EPANET .inp file.

    Returns
    -------
    G : nx.Graph
        NetworkX graph of the water network.
    """
    wn = epanet(inp_path)

    # Node metadata
    node_ids = wn.getNodeNameID()
    node_types = wn.getNodeType()
    node_demands = next(iter(wn.getNodeBaseDemands().values()))
    node_elevs = wn.getNodeElevations()
    node_coords = list(zip(wn.getNodeCoordinates()['x'].values(), wn.getNodeCoordinates()['y'].values()))
    G = nx.Graph()

    for i, node_id in enumerate(node_ids):
        node_type = node_types[i]
        demand = float(node_demands[i])
        elev = float(node_elevs[i])
        x, y = node_coords[i]

        is_source = node_type in ("RESERVOIR", "TANK")
        is_consumer = (node_type == "JUNCTION" and demand > 0)

        G.add_node(
            i+1,
            node_id=node_id,
            type=node_type,
            demand=demand,
            elevation=elev,
            pos=np.array([x, y]),
            is_source=is_source,
            is_consumer=is_consumer,
            weight=demand if demand > 0 else 0.0,
        )

    # Edge metadata
    link_ids = wn.getLinkNameID()
    link_types = wn.getLinkType()
    link_nodes = wn.getLinkNodesIndex()
    link_len = wn.getLinkLength()
    link_diam = wn.getLinkDiameter()
    link_rough = wn.getLinkRoughnessCoeff()
    # link_type_map = {0: "pipe", 1: "pump", 2: "valve"}

    for i, link_id in enumerate(link_ids):
        n1, n2 = link_nodes[i]
        ltype = link_types[i]

        length = float(link_len[i]) if i < len(link_len) else None
        diam = float(link_diam[i]) if i < len(link_diam) else None
        rough = float(link_rough[i]) if i < len(link_rough) else None

        G.add_edge(
            n1,
            n2,
            id=link_id,
            type=ltype,
            length=length,
            diameter=diam,
            roughness=rough,
            weight=length if length else None,
        )

    wn.closeNetwork()
    G = nx.relabel_nodes(G, dict(zip(G.nodes, range(len(G)))))
    G.graph["sources"] = [node for node, data in G.nodes.items() if data.get("is_source")]
    return G

def get_all_water_networks(names: list=None) -> Dict[str, nx.Graph]:
    """
    Scan subfolders of data/water for .inp files (excluding names containing '_temp'),
    convert each found .inp to a NetworkX graph using water_inp_to_nx and return a dict
    mapping folder name -> graph.
    """
    networks = {}
    base_folder = r"data\water"

    if not os.path.isdir(base_folder):
        return networks

    networks_list = names if names else os.listdir(base_folder)

    for folder_name in networks_list:
        folder_path = os.path.join(base_folder, folder_name)

        if not os.path.isdir(folder_path):
            continue

        inp_files = [
            fn for fn in os.listdir(folder_path)
            if fn.lower().endswith(".inp") and "_temp" not in fn
        ]

        if not inp_files:
            continue

        inp_path = os.path.join(folder_path, inp_files[0])
        try:
            networks[folder_name] = water_inp_to_nx(inp_path)
        except Exception as e:
            # keep behavior simple: skip on error
            print(f"Warning: failed to load {inp_path}: {e}")

    return networks

############################
#          power           #
############################
    

def get_sfo_networks(improve: bool=False) -> Dict[str, nx.Graph]:
    networks_dict = {
        'SFO Davidson': get_sfo_bounded_network("P3U", "box"),
        'SFO Pacific': get_sfo_bounded_network("P2U", "box2"),
    }
    if improve:        
        G = networks_dict['SFO Pacific']
        networks_dict['SFO Pacific'], _ = improve_tree(G, R=8, root=G.graph["sources"][0], weight_quantile=0.3, show_progress=True)

        G = networks_dict['SFO Davidson']
        networks_dict['SFO Davidson'], _ = improve_tree(G, R=6, root=G.graph["sources"][0], weight_quantile=0.3, show_progress=True)
    return networks_dict

def get_sfo_bounded_network(area: str="P2U", box: str="box", include_ties: bool=True) -> nx.Graph:
    dir = fr"data\power\better_grids\SFO\{area}"
    boundery = gpd.read_file(fr"{dir}\{box}.geojson").to_crs(3857)
    G = build_mv_network_from_boundery(dir, boundery.to_crs(32610), full_network=False)
    G = _radial_graph_with_ties(G, include_ties=include_ties)
    _change_graph_crs(G, old_crs=32610, new_crs=3857)
    G = nx.relabel_nodes(G, dict(zip(G.nodes, range(len(G)))))
    G.graph["sources"] = [node for node, data in G.nodes.items() if data["is_source"]]

    return G


def _radial_graph_with_ties(G: nx.Graph, include_ties: bool=True) -> nx.Graph:
    normally_closed_edges = [
        (u, v) for u, v, data in G.edges(data=True)
        if not data.get("is_tie", False)
    ]
    closed_graph = nx.edge_subgraph(G, normally_closed_edges).copy()
    radial_edges = list(nx.minimum_spanning_tree(closed_graph, weight="length_m").edges)
    selected_edges = set(map(tuple, radial_edges))
    if include_ties:
        selected_edges.update((u, v) for u, v, data in G.edges(data=True) if data.get("is_tie", False))
    return nx.edge_subgraph(G, selected_edges).copy()


def _change_graph_crs(G: nx.Graph, old_crs: int, new_crs: int):
    """inplace change the pos attribute of the graph"""
    pos = nx.get_node_attributes(G, "pos").items()
    s = gpd.points_from_xy(list(map(lambda x: x[1][0], pos)), list(map(lambda x: x[1][1], pos)), crs=old_crs).to_crs(new_crs)
    for (node, _), pos in zip(pos, s):
        G.nodes[node]["pos"] = (pos.x, pos.y)

def build_network_from_gdfs(
    lines: gpd.GeoDataFrame,
    transf: gpd.GeoDataFrame,
    subs: gpd.GeoDataFrame,
):
    """Construct a networkx Graph from GeoDataFrames of lines, transformers and substations.

    Expects `lines` to contain LineString geometries and columns NodeA/NodeB
    identifying endpoints. Transformer and substation GeoDataFrames should
    contain Node identifiers and point geometries.
    """
    # --- Prepare graph ---
    G = nx.Graph()

    # --- Add line edges ---
    for _, row in lines.iterrows():
        geom = row.geometry
        if not isinstance(geom, LineString):
            continue
        # Edge endpoints (as coordinates)
        start, end = geom.coords[0], geom.coords[-1]
        nodeA, nodeB = _clean_node_id(row.NodeA), _clean_node_id(row.NodeB)
        length = geom.length  # Euclidean length in same units as CRS
        edge_attrs = _line_row_edge_attrs(row, length)
        G.add_edge(nodeA, nodeB, geometry=geom, **edge_attrs)
        G.nodes[nodeA]["pos"] = np.array(list(start))
        G.nodes[nodeB]["pos"] = np.array(list(end))
        is_switch_boundary = edge_attrs["is_switch"] or edge_attrs["is_tie"]
        G.nodes[nodeA]["middle_node"] = G.nodes[nodeA].get("middle_node", True) and not is_switch_boundary
        G.nodes[nodeB]["middle_node"] = G.nodes[nodeB].get("middle_node", True) and not is_switch_boundary
        for node in (nodeA, nodeB):
            G.nodes[node].setdefault("is_source", False)
            G.nodes[node].setdefault("weight_kva", 0)

    # --- Add transformers as load nodes ---
    for _, row in transf.iterrows():
        pt = row.geometry
        node = _clean_node_id(row.Node)
        G.add_node(node)
        G.nodes[node]["weight_kva"] = row.get("Size_KVA", 0)
        G.nodes[node]["is_source"] = False
        G.nodes[node]["pos"] = np.array([row.geometry.x, row.geometry.y])
        G.nodes[node]["middle_node"] = False

    # --- Mark substations as sources ---
    for _, row in subs.iterrows():
        pt = row.geometry
        pos = np.array([pt.x, pt.y])
        node = _clean_node_id(row.Node)
        if node not in G:
            G.add_node(node)
        G.nodes[node]["is_source"] = True
        G.nodes[node]["weight_kva"] = 0
        G.nodes[node]["pos"] = pos
        G.nodes[node]["middle_node"] = False

    # --- remove middle nodes ---
    G = _remove_middle_nodes(G)
    G = nx.subgraph(G, max(nx.connected_components(G), key=len))
    G.graph["sources"] = set([node for node, data in G.nodes.items() if data["is_source"]])
    return G


def _clean_node_id(node) -> str:
    return str(node).strip()


def _edge_key_from_nodes(node_a, node_b) -> tuple[str, str]:
    return tuple(sorted((_clean_node_id(node_a), _clean_node_id(node_b))))


def _line_row_edge_attrs(row, length: float) -> dict:
    is_tie = bool(getattr(row, "is_tie", False))
    is_switch = bool(getattr(row, "is_switch", False)) or is_tie
    status = getattr(row, "Status", None)
    normally_closed = bool(getattr(row, "normally_closed", not is_tie))
    return {
        "length_m": float(length),
        "length": float(length),
        "is_switch": is_switch,
        "is_tie": is_tie,
        "normally_closed": normally_closed,
        "status": status,
        "code": getattr(row, "Code", None),
        "switch_kind": getattr(row, "switch_kind", None),
        "switch_codes": getattr(row, "switch_codes", None),
    }


def _add_switch_metadata_to_lines(
    lines: gpd.GeoDataFrame,
    switches: gpd.GeoDataFrame | None = None,
) -> gpd.GeoDataFrame:
    """Annotate SMART-DS lines with switch and normally-open tie metadata."""
    lines = lines.copy()
    lines["NodeA"] = lines["NodeA"].map(_clean_node_id)
    lines["NodeB"] = lines["NodeB"].map(_clean_node_id)
    lines["_edge_key"] = lines.apply(lambda row: _edge_key_from_nodes(row.NodeA, row.NodeB), axis=1)

    if "Status" in lines:
        status = lines["Status"].astype(str).str.strip()
        open_edges = set(lines.loc[status == "0", "_edge_key"])
    else:
        open_edges = set()

    switch_edges = set()
    switch_kind_by_edge = {}
    switch_codes_by_edge = {}
    if switches is not None and len(switches) > 0:
        switches = switches.copy()
        switches["NodeA"] = switches["NodeA"].map(_clean_node_id)
        switches["NodeB"] = switches["NodeB"].map(_clean_node_id)
        switches["_edge_key"] = switches.apply(lambda row: _edge_key_from_nodes(row.NodeA, row.NodeB), axis=1)
        switches["switch_kind"] = switches["Code"].astype(str).str.extract(r"^([^\(]+)")[0].str.strip()
        switch_edges = set(switches["_edge_key"])
        switch_kind_by_edge = switches.groupby("_edge_key")["switch_kind"].apply(
            lambda values: "|".join(sorted(set(map(str, values))))
        ).to_dict()
        switch_codes_by_edge = switches.groupby("_edge_key")["Code"].apply(
            lambda values: "|".join(sorted(set(map(str, values))))
        ).to_dict()

    lines["is_switch"] = lines["_edge_key"].isin(switch_edges | open_edges)
    lines["is_tie"] = lines["_edge_key"].isin(open_edges)
    lines["normally_closed"] = ~lines["is_tie"]
    lines["switch_kind"] = lines["_edge_key"].map(switch_kind_by_edge)
    lines["switch_codes"] = lines["_edge_key"].map(switch_codes_by_edge)
    return lines.drop(columns=["_edge_key"])


def build_mv_network(data_dir, voltage=12.47):
    """
    Build a NetworkX graph from SMART-DS shapefiles for a given voltage level (e.g., 12.47 kV).
    
    Parameters
    ----------
    data_dir : str
        Path to the folder containing Line_N, DistribTransf_N, and HVMVSubstation_N shapefiles.
    voltage : float
        Nominal voltage (default: 12.47).
    
    Returns
    -------
    G : nx.Graph
        Network graph with attributes:
        - edges: geometry, length_m
        - nodes: weight (Size_KVA), is_source
    """

    # --- Load shapefiles ---
    lines = gpd.read_file(f"{data_dir}/Line_N.shp")
    transf = gpd.read_file(f"{data_dir}/DistribTransf_N.shp")
    subs = gpd.read_file(f"{data_dir}/HVMVSubstation_N.shp")
    switches = _read_switching_devices(data_dir)

    # --- Filter lines by voltage ---
    vcol = [c for c in lines.columns if "NomV" in c or "voltage" in c][0]
    lines = lines[np.isclose(lines[vcol].astype(float), voltage, atol=0.1)]
    if switches is not None and "NomV_kV" in switches:
        switches = switches[np.isclose(switches["NomV_kV"].astype(float), voltage, atol=0.1)]
    lines = _add_switch_metadata_to_lines(lines, switches)

    return build_network_from_gdfs(
        lines, transf, subs
    )



def build_mv_network_from_boundery(data_dir, boundery, full_network: bool=True):
    """
        Build a network of sum sub region in SFO. Takes the lines as lines,
        transformere as targets and the substations as sources.
        load only the data in the boundery polygon
        If full_network == True, the lines are all the optinal lines, and other wise, its the real lines
    """
    # Read the shapefile
    if full_network:
        lines = gpd.read_file(f"{data_dir}/StreetMap_branches.shp").set_crs(32610)
        lines.rename(columns={'Node_A': 'NodeA', 'Node_B': 'NodeB'}, inplace=True)
        nodes = gpd.read_file(f"{data_dir}/StreetMap_nodes.shp").set_crs(32610)
        nodes.rename(columns={'Node_A': 'NodeA', 'Node_B': 'NodeB'}, inplace=True)
    else:
        lines = gpd.read_file(f"{data_dir}/Line_N.shp").set_crs(32610)
        nodes = gpd.read_file(f"{data_dir}/DistribTransf_N.shp").set_crs(32610)
    switches = _read_switching_devices(data_dir)
    if switches is not None:
        switches = switches.set_crs(32610)
    transf = gpd.read_file(f"{data_dir}/DistribTransf_N.shp").set_crs(32610)
    subs = gpd.read_file(f"{data_dir}/HVMVSubstation_N.shp").set_crs(32610)
    
    # filter by polygon
    lines = lines.sjoin(boundery, how='inner')
    if switches is not None:
        switches = switches.sjoin(boundery, how='inner')
    transf = transf.sjoin(boundery, how='inner')
    subs = subs.sjoin(boundery, how='inner')
    nodes = nodes.sjoin(boundery, how='inner')
    lines = _add_switch_metadata_to_lines(lines, switches)
    
    # get closest transf_for_each_node
    # ensure consistent integer index positions
    transf = _find_closest_node2_to_node1(transf, nodes)
    subs = _find_closest_node2_to_node1(subs, nodes)

    return build_network_from_gdfs(
        lines, transf, subs
    )


def _read_switching_devices(data_dir: str) -> gpd.GeoDataFrame | None:
    switch_path = os.path.join(data_dir, "SwitchingDevices_N.shp")
    if not os.path.exists(switch_path):
        return None
    return gpd.read_file(switch_path)

def _find_closest_node2_to_node1(nodes1: gpd.GeoDataFrame, nodes2: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Map each row in nodes1 to the closest row in nodes2 by geometry.

    Returns a copy of nodes1 with the 'Node' column replaced by the closest
    Node value from nodes2. Uses a KD-tree for efficiency.
    """
    nodes2 = nodes2.reset_index(drop=True)
    nodes1 = nodes1.reset_index(drop=True)

    # extract point coordinates (use centroid for safety)
    nodes2_coords = np.column_stack((nodes2.geometry.x.values, nodes2.geometry.y.values))
    nodes1_coords = np.column_stack((nodes1.geometry.x.values, nodes1.geometry.y.values))

    # build kd-tree and query nearest node for each transformer
    tree = cKDTree(nodes2_coords)
    _, idx = tree.query(nodes1_coords, k=1)

    # assign transformer's Node to the closest node's Node value
    nodes1["Node"] = nodes2.loc[idx, "Node"].values
    return nodes1



def _get_boundary_nodes(G: nx.Graph, subset: set) -> set:
    """Return nodes adjacent to `subset` but not in `subset`."""
    neighbors = set()
    for node in subset:
        neighbors.update(G.neighbors(node))
    return neighbors - subset


def _remove_middle_nodes(G: nx.Graph) -> nx.Graph:
    """
    Remove nodes marked as 'middle_node=True' and reconnect their boundary nodes.

    Parameters
    ----------
    G : nx.Graph
        Original graph with node attribute 'middle_node' (bool).

    Returns
    -------
    G_clean : nx.Graph
        Simplified graph without middle nodes, where
        boundary nodes are directly connected.
    """
    # Identify all middle nodes
    middle_nodes = {n for n, d in G.nodes(data=True) if d.get("middle_node", False)}

    # Build subgraph of middle nodes
    G_middle = G.subgraph(middle_nodes)

    # Start with subgraph of non-middle nodes
    G_clean = G.subgraph(set(G.nodes) - middle_nodes).copy()

    # Reconnect boundaries across removed middle nodes
    for comp in nx.connected_components(G_middle):
        boundary = list(_get_boundary_nodes(G, comp))
        if len(boundary) < 2:
            continue
        comp_with_boundary = set(comp) | set(boundary)
        local_graph = G.subgraph(comp_with_boundary)
        for i, node in enumerate(boundary):
            for b in boundary[i + 1:]:
                if node == b:
                    continue
                try:
                    path = nx.shortest_path(local_graph, node, b, weight="length_m")
                except nx.NetworkXNoPath:
                    continue
                G_clean.add_edge(node, b, **_aggregate_path_attrs(G, path))

    return G_clean


def _aggregate_path_attrs(G: nx.Graph, path: list) -> dict:
    edges = list(zip(path[:-1], path[1:]))
    attrs = [G.edges[e] for e in edges]
    is_tie = any(bool(data.get("is_tie", False)) for data in attrs)
    is_switch = any(bool(data.get("is_switch", False)) for data in attrs)
    length_m = sum(float(data.get("length_m", data.get("length", 0.0))) for data in attrs)
    switch_kind = sorted({
        kind
        for data in attrs
        for kind in str(data.get("switch_kind", "")).split("|")
        if kind and kind != "nan" and kind != "None"
    })
    switch_codes = sorted({
        code
        for data in attrs
        for code in str(data.get("switch_codes", "")).split("|")
        if code and code != "nan" and code != "None"
    })
    return {
        "length_m": length_m,
        "length": length_m,
        "is_switch": is_switch,
        "is_tie": is_tie,
        "normally_closed": not is_tie,
        "status": 0 if is_tie else 1,
        "switch_kind": "|".join(switch_kind) if switch_kind else None,
        "switch_codes": "|".join(switch_codes) if switch_codes else None,
        "original_edges": [tuple(edge) for edge in edges],
        "original_nodes": list(path),
    }




def get_all_power_networks() -> Dict[str, nx.Graph]:
    networks_names = ["iceland"]
    # make dict
    power_net_dict = {}
    for name in networks_names:
        print(name)
        try:
            graph = from_pp(name)
            power_net_dict[name] = graph
        except:
            continue
    # make df
    power_df = networks_dict_to_df(power_net_dict)
    return power_df


def from_pp(network_name) -> nx.Graph:
    network = eval(f"pp.networks.{network_name}()")

    G = pp.topology.create_nxgraph(network)
    pos = {
        node: np.array(list(shape(json.loads(node_pos)).centroid.coords)).reshape(-1)
        for node, node_pos in network.bus.geo.items()
    }
    nx.set_node_attributes(G, pos, "pos")
    G.graph["sources"] = set(network.gen.bus) | set(network.ext_grid.bus)
    nx.set_node_attributes(G, {node: node in G.graph["sources"] for node in G.nodes}, "is_source")
    return G




def read_dss(folder: str):
    # --- Helper to parse key=value pairs in DSS lines ---
    def parse_params(line):
        params = {}
        for match in re.finditer(r'(\w+)=([\w\.\-\+]+)', line):
            key, val = match.groups()
            params[key.lower()] = val
        return params

    # --- 1. Read Bus coordinates ---
    bus_pos = {}
    with open(f"{folder}/Buscoords.dss") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                bus, x, y = parts[:3]
                # print(bus, x, y)
                bus_pos[bus.lower()] = (float(x), float(y))

    # --- 2. Read Lines ---
    edges = []
    with open(f"{folder}/Lines.dss") as f:
        for line in f:
            if not line.strip().lower().startswith("new line"):
                continue
            p = parse_params(line)
            b1 = p.get("bus1", "").split(".")[0].lower()
            b2 = p.get("bus2", "").split(".")[0].lower()
            if b1 and b2:
                edges.append((b1, b2, p))

    # --- 3. Read Loads ---
    loads = {}
    with open(f"{folder}/Loads.dss") as f:
        for line in f:
            if not line.strip().lower().startswith("new load"):
                continue
            p = parse_params(line)
            bus = p.get("bus1", "").split(".")[0].lower()
            kw = float(p["kw"]) if "kw" in p else None
            loads[bus] = kw

    # --- 4. Find source bus ---
    source_bus = None
    with open(f"{folder}/Master.dss") as f:
        for line in f:
            if line.lower().startswith("new circuit"):
                m = re.search(r"bus1\s*=\s*([\w\-\_]+)", line.lower())
                if m:
                    source_bus = m.group(1).lower()
                break

    # --- 5. Build the NetworkX graph ---
    G = nx.Graph()

    # add edges
    for b1, b2, attr in edges:
        G.add_edge(b1, b2, **attr)

    # add node attributes: position, load, source flag
    for n in list(G.nodes):
        if bus_pos.get(n) is None:
            G.remove_node(n)
        else:
            G.nodes[n]["pos"] = np.array(bus_pos.get(n))
            G.nodes[n]["load_kw"] = loads.get(n)
            G.nodes[n]["is_source"] = (n == source_bus)

    print(nx.info(G))
    print("Example node:", list(G.nodes(data=True))[:3])

    draw_network(G, with_labels=False, strc_size=10, chain_size=0)

    def get_max2copms(G):
        comps2 = [x for comp in nx.k_edge_components(G, 2) if len(comp) > 1 for x in comp]
        return nx.subgraph(G, comps2).copy()

    # plt.figure()
    draw_network(get_max2copms(G), with_labels=False, strc_size=10, chain_size=0, width=4)


##### SFO ####
import boto3
from botocore import UNSIGNED
from botocore.config import Config
import os
from scipy.spatial import cKDTree

def download_sfo_data(area: str):
    """ Download data from better grids repo"""
    # 1️⃣ Configure S3 public client
    s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))
    area = 'P3R'
    bucket = "oedi-data-lake"
    prefix = f"SMART-DS/v1.0/GIS/SFO/{area}/"

    # 2️⃣ Where to save locally
    local_dir = rf"data\power\better_grids\SFO\{area}"
    os.makedirs(local_dir, exist_ok=True)

    # 3️⃣ List and download all files under prefix
    paginator = s3.get_paginator("list_objects_v2")
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        for obj in page.get("Contents", []):
            key = obj["Key"]
            rel_path = key[len(prefix):]  # relative path inside the folder
            local_path = os.path.join(local_dir, rel_path)

            # Create any missing folders
            os.makedirs(os.path.dirname(local_path), exist_ok=True)

            # Skip empty "folders"
            if key.endswith("/"):
                continue

            print(f"⬇️ Downloading {key} → {local_path}")
            s3.download_file(bucket, key, local_path)

    print("\n✅ Done! All files saved in:", local_dir)
