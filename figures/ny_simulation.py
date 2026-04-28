import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from shapely.geometry import Polygon, Point, box, MultiPolygon
import geopandas as gpd
import indexes.simulation as sm
from scipy.spatial.distance import squareform, pdist
from scipy.spatial import Voronoi
from utilities.lkh_tsp import solve_tsp_lin_kernighan
from indexes import GraphRel, edge_probs_by_length, Float
from typing import Literal
import utilities.read_write as rw
from utilities.figures_utilities import saidi_with_lengths, failing_rate_from_spanning_tree
from figures.ny_plot import NEW_YORK_FILE_NAME, NEW_YORK_ASYMPTOTIC_FILE_NAME, NEW_YORK_MCMC_FILE_NAME
from indexes.simulation import simulate_rel



#### The main functions to generate data for the graph ####
def ny_graph_simulation(file_name: str=NEW_YORK_FILE_NAME) -> pd.DataFrame:
    """
    This is the main function used to generate the simulation for the
    ny graph. Here we take one network with 600 nodes, start from a base cycle
    and add random redundant edges
    the SAIDI is computed using MCMC Simulations
    """
    manhattan_polygon = get_polygon_for_simulations()
    p_list = sorted(
        [0.0, 5e-6, 7e-6] +
        list(np.round(np.arange(2e-5, 1e-4+5e-6, 1e-5), 6)) +
        list(np.round(np.arange(1e-5, 1e-3+5e-5, 5e-5), 5))
    )
    print(f"{len(p_list)=}")

    # Primary RNG for the entire run — not optional/None. Change the integer
    # here to get a different deterministic outcome.
    rng = np.random.default_rng(4000)

    df = generate_and_improve_network_by_polygon(
        polygon=manhattan_polygon,
        n_trials=20,
        n_points=600,
        n_sources=1,
        r_list=list(range(0, 21)) + [301],
        p_list=p_list,
        rng=rng,
        add_tree=True,
        file_name=file_name,
        time_limit=15,
    )
    return df

def ny_asymptotic_graph_simulation(file_name: str=NEW_YORK_ASYMPTOTIC_FILE_NAME) -> pd.DataFrame:
    """Generate asymptotic simulations by increasing network size per slice.

    Uses the same deterministic RNG pattern as ``ny_graph_simulation``.
    """
    nyc__polygon = get_polygon_for_simulations()
    p_list = [5e-4]

    # use the same style: define a top-level RNG and pass it downward
    rng = np.random.default_rng(5000)

    df = generate_and_improve_networks_by_polygon(
        polygon=nyc__polygon,
        n_trials=5,
        n_networks=10,
        n_points_per_iter=[100 + 30 * i ** 2 for i in range(0, 10)],
        n_points_per_source=600,
        r_list_per_iter=[[i + 3] for i in range(0, 10)],
        p_list=p_list,
        add_tree=False,
        file_name=file_name,
        rng=rng,
        time_limit=20,
    )
    return df



def get_polygon_for_simulations() -> Polygon:
    """
    Load Manhattan polygon from GeoJSON files.
    
    Returns
    -------
    Polygon
        Manhattan polygon
    """
    # Load Manhattan polygon from GeoJSON
    manhattan_poly_file = r"data\ny\manhattan.geojson"
    gdf = gpd.read_file(manhattan_poly_file).set_crs(4326)
    gdf['area'] = gdf.area
    manhattan_polygon = gdf.to_crs(4326).loc[gdf.area.idxmax()].geometry
    # nyc
    return manhattan_polygon


##### MCMC simulation ####
def generate_data_MCMC(
    file_name: str=NEW_YORK_MCMC_FILE_NAME,
    manhattan_file_name: str=NEW_YORK_FILE_NAME,
)->pd.DataFrame:
    """
    Generate reliability simulation data for plotting.
    
    Loads a network from the Manhattan dataset and simulates reliability
    for multiple p values, then saves the results to a CSV file.
    """    
    ny_data = rw.read_nxjson(manhattan_file_name)
    res_data = [] # r, p, t, saidi
    rng = np.random.default_rng(100)
    for r in [0, 1]:
        graph, sources = ny_data.query('r == @r').iloc[0][["graph", "sources"]].values
        for p_mean in [5e-4, 1e-3, 5e-3]:
            # set edge probs
            edges_prob, _ = edge_probs_by_length(graph, p=p_mean, mode="mean")
            edges_prob = {e: Float(p) for e, p in edges_prob.items()}

            # set GraphRel to handle multi sources
            gr = GraphRel(graph, edges_prob=edges_prob, sources=sources)

            # simulate
            sim_res, _ = simulate_rel(
                nx.Graph(gr.graph),
                source=gr.source,
                prob_attr="prob",
                rel_type="saidi",
                weight_attr="weight",
                T_days=365,
                mean_cycle_days=0.5,
                rng=rng,
                show_progress=True
            )
            # combine results to data
            times = sim_res.times
            cum_saidis = sim_res.cum_rel
            p_means = [p_mean] * len(times) 
            rs = [r] * len(times) 
            current_res_data = list(zip(rs, p_means, times, cum_saidis))
            res_data.extend(current_res_data)
    res_df: pd.DataFrame = pd.DataFrame(res_data, columns=["r", "p", "t", "saidi"])
    rw.write_nxjson(res_df, file_name)
    print("DONE!")
    return res_df



####  generate networks for polygon###
def generate_and_improve_networks_by_polygon(
    polygon: Polygon,
    n_trials: int,
    n_networks: int,
    n_points_per_iter: list[int],
    n_points_per_source: int,
    r_list_per_iter: list[list[int]],
    p_list: list[float],
    add_tree: bool,
    file_name: str = None,
    rng: np.random.Generator = None,
    time_limit: float = None,
) -> pd.DataFrame:
    """
    Generate and improve networks by slicing a polygon into horizontal tiers.
    
    Subdivides a geographic polygon into n_networks horizontal slices, then
    generates and improves networks for each slice with potentially different
    node counts and redundancy targets per tier.
    
    Parameters
    ----------
    polygon : Polygon
        Geographic region to subdivide and generate networks in.
    n_trials : int
        Number of trials per slice/tier.
    n_networks : int
        Number of horizontal slices to create.
    n_points_per_iter : list[int]
        Number of nodes for each tier (must have length n_networks).
    n_points_per_source : int
        Divisor for computing number of sources (n_sources = ceil(n_points / n_points_per_source)).
    r_list_per_iter : list[int]
        Target redundancy values for each tier (must have length n_networks).
    p_list : list[float]
        Failure probabilities for reliability evaluation.
    add_tree : bool
        Include MST variants in addition to rings and improved networks.
    file_name : str, optional
        Output JSON file path. If provided, saves results.
    time_limit : float, optional
        TSP solver time limit in seconds.
    
    Returns
    -------
    pd.DataFrame
        Combined results from all tiers.
        Columns: n, m, r, rho, sources, p, rel, graph, trial, iter.
    
    Raises
    ------
    ValueError
        If n_points_per_iter or r_list_per_iter lengths don't match n_networks.
    """
    # validate
    if len(n_points_per_iter) != n_networks:
        raise ValueError(f"Length of n_points_list {len(n_points_per_iter)} should match the number of networks {n_networks}")
    if len(r_list_per_iter) != n_networks:
        raise ValueError(f"Length of r_list {len(r_list_per_iter)} should match the number of networks {n_networks}")
    # prepare polygon
    valid_polygon = polygon.buffer(0)
    min_x, min_y, max_x, max_y = valid_polygon.bounds
    it = 0
    all_data = []
    for y in np.linspace(min_y, max_y, n_networks + 1)[:-1][::-1]:
        print(f'----- {it} -----')
        # slice polygon
        bbox = box(min_x, y, max_x, max_y)
        sliced_polygon = valid_polygon.intersection(bbox)
        # generate network for the sliced polygon
        n_points = n_points_per_iter[it]
        n_sources = max([1, int(np.ceil(n_points / n_points_per_source))])
        # create a per-iteration RNG derived from the top-level RNG for reproducibility
        iter_rng = np.random.default_rng(rng.integers(0, 2 ** 31 - 1)) if rng is not None else np.random.default_rng()
        iter_data = generate_and_improve_network_by_polygon(
            polygon=sliced_polygon,
            n_trials=n_trials,
            n_points=n_points,
            n_sources=n_sources,
            r_list=r_list_per_iter[it],
            p_list=p_list,
            add_tree=add_tree,
            rng=iter_rng,
            time_limit=time_limit,
        )
        # add to all data
        iter_data["iter"] = it
        all_data.append(iter_data)
        it += 1
    # save data
    untied_data = pd.concat(all_data)
    if file_name is not None:
        rw.write_nxjson(untied_data, file_name)
        print(f"save data to {file_name}")
    return untied_data


def generate_and_improve_network_by_polygon(
    polygon: Polygon,
    n_trials: int,
    n_points: int,
    n_sources: int,
    r_list: list[int],
    p_list: list[float],
    add_tree: bool,
    file_name: str = None,
    rng: np.random.Generator = None,
    time_limit: float = None,
) -> pd.DataFrame:
    """
    Generate and improve networks for a single polygon with multiple trials.
    
    Creates base ring networks, optionally adds MST variants, applies edge
    improvements to reach target redundancies, computes reliability metrics,
    and optionally saves results to a JSON file.
    
    Parameters
    ----------
    polygon : Polygon
        Geographic region to generate networks in.
    n_trials : int
        Number of independent trials/base rings to generate.
    n_points : int
        Number of nodes per network.
    n_sources : int
        Number of source/supply nodes.
    r_list : list[int]
        Target redundancy values to save improved graphs for.
    p_list : list[float]
        Failure probability values for reliability evaluation.
    add_tree : bool
        Whether to include MST variant in addition to ring and improved nets.
    file_name : str, optional
        Output JSON file path. If provided, saves results.
    rng: np.random.Generator
        Random generator
    time_limit : float, optional
        TSP solver time limit in seconds.
    
    Returns
    -------
    pd.DataFrame
        One row per (trial, network, p_value) combination.
        Columns: n, m, r, rho, sources, p, rel, graph, trial.
    """
    graphs_data = pd.DataFrame()
    for trial in range(n_trials):
        # derive a per-trial RNG from the top-level RNG to keep runs reproducible
        trial_rng = np.random.default_rng(rng.integers(0, 2 ** 31 - 1)) if rng is not None else np.random.default_rng()
        networks: list[nx.Graph] = []
        # base ring
        base_ring = base_ring_from_polygon(
            polygon=polygon,
            num_points=n_points,
            n_sources=n_sources,
            rng=trial_rng,
            time_limit=time_limit,
        )
        if 1 in r_list:
            networks.append(base_ring)

        # add tree
        if add_tree:
            points = list(nx.get_node_attributes(base_ring, "pos").values())
            tree = _minimum_spanning_tree(points, rng=trial_rng)
            tree.graph["sources"] = base_ring.graph["sources"]
            networks.append(tree)

        # add improved nets
        improved_networks = random_improve_network(
            base_graph=base_ring,
            r_list=r_list,
            rng=trial_rng,
        )
        networks.extend(improved_networks.values())

        # tranform to df with attributes
        graphs_data_temp = _graphs_df_from_graphs_list(networks, p_list=p_list, rng=trial_rng)
        graphs_data_temp['trial'] = trial
        graphs_data = pd.concat([graphs_data, graphs_data_temp])
    if file_name:
        rw.write_nxjson(graphs_data, file_name)
    return graphs_data


# Calculate network attributes
def _graphs_df_from_graphs_list(graphs_list: list[nx.Graph], p_list: list[float], rng: np.random.Generator = None) -> pd.DataFrame:
    """
    Convert a list of graphs into a DataFrame with computed network metrics.
    
    For each graph, computes: n (nodes), m (edges), r (redundancy), rho (r/n),
    sources list, and reliability/SAIDI values across all p values.
    
    Parameters
    ----------
    graphs_list : list[nx.Graph]
        List of NetworkX graphs to process.
    p_list : list[float]
        Failure probability values.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: n, m, r, rho, sources, p, rel, graph.
        Sorted by redundancy (r).
    
    Notes
    -----
    Graphs with rho > 0.4 are skipped (set rel to 0) as they exceed
    practical network density assumptions.
    """
    graphs_data = pd.DataFrame(
        columns=["n", "m", "r", 'rho', "sources", "p", "p_rate", "rel", "graph"]
    )
    for graph in tqdm(graphs_list, desc="Calc SAIDI"):
        n = len(graph)
        m = len(graph.edges)
        r = m - n + 1
        rho = r / n
        sources = graph.graph["sources"]
        # derive a graph-specific RNG for the simulation step to avoid
        # exhausting the trial RNG and to keep deterministic reproducibility
        graph_rng = np.random.default_rng(rng.integers(0, 2 ** 31 - 1)) if rng is not None else np.random.default_rng()
        saidi_rate_by_p = _saidi_using_simulation(graph, p_list, sources, rng=graph_rng)
        saidis = [saidi_rate_by_p[p]["saidi"] for p in p_list]
        p_rates = [saidi_rate_by_p[p]["p_rate"] for p in p_list]
        graphs_data.loc[len(graphs_data)] = {
            'n': n, 'm': m, 'r': r, 'rho': rho,     
            'sources': sources, 'p': list(p_list), 'p_rate': p_rates,
            'rel': list(saidis), 'graph': graph
        }
    graphs_data.sort_values('r', inplace=True)
    return graphs_data

def _saidi_using_simulation(
        G: nx.Graph,
        p_list: list[float],
        sources: list,
        rng: np.random.Generator,
        p_rate: float | None = None,
    ) -> dict[float, dict[Literal["p_rate", "saidi"], float]]:
    """
    Compute SAIDI values for a list of failure probabilities.
    
    Parameters
    ----------
    G : nx.Graph
        The network graph.
    p_list : list[float]
        List of failure probabilities.
    sources : list
        List of source nodes.
    rng : np.random.Generator
        Random number generator.
    p_rate : float | None, optional
        Pre-computed p_rate.
        
    Returns
    -------
    dict
        Dictionary mapping p to computed p_rate and SAIDI values.
    """
    res_dict = {}
    for p in p_list:
        # infer the p_rate from the input mean_p
        p_rate = edge_probs_by_length(G, p=p, mode="mean", tol=1e-6, max_iter=50)[1]
        saidi = saidi_with_lengths(G, sources=sources, p=p_rate, mode="rate", rng=rng, show_progress=False, mean_cycle_days=0.1)
        res_dict[p] = {"p_rate": p_rate, "saidi": saidi}
    return res_dict


##### improve network #####

def _get_random_candidate_edges(G: nx.Graph, rng: np.random.Generator) -> list[tuple]:
    """Return a randomly shuffled list of candidate edges using `rng`.

    Parameters
    ----------
    G : nx.Graph
        Input graph
    rng : np.random.Generator
        Required RNG to draw permutations from.

    Returns
    -------
    list[tuple]
        Randomly shuffled list of candidate edges
    """
    nodes = list(G.nodes)

    # Get all possible edges that don't exist in the graph
    candidate_edges = [(u, v) for i, u in enumerate(nodes) for j, v in enumerate(nodes) if i < j and not G.has_edge(u, v)]

    if len(candidate_edges) == 0:
        return []

    # Shuffle randomly using provided generator
    shuffled_edges = list(rng.permutation(candidate_edges))

    return shuffled_edges


def random_improve_network(
        base_graph: nx.Graph,
        r_list: list[int],
        rng: np.random.Generator,
    ) -> dict:
    """Improve network by adding random edges drawn using `rng`.

    Parameters
    ----------
    base_graph : nx.Graph
        Base network graph
    r_list : list[int]
        List of target redundancy values to save graphs for
    rng : np.random.Generator
        Required RNG for random choices and shuffling

    Returns
    -------
    dict
        Dictionary mapping redundancy r to improved graph copies
    """
    # Get candidate edges (randomly shuffled by provided generator)
    candidate_edges = _get_random_candidate_edges(base_graph, rng=rng)

    if len(candidate_edges) == 0:
        return {}

    # init vars
    graph = base_graph.copy()
    max_r = max(r_list)
    r = len(graph.edges) - len(graph.nodes) + 1
    improved_graphs = {}

    # Randomly select from candidate edges one by one
    candidate_idx = 0

    # add random edges
    while r < max_r and candidate_idx < len(candidate_edges):
        u, v = candidate_edges[candidate_idx]
        candidate_idx += 1

        # check if non existing and connect 2-degree nodes
        if not graph.has_edge(u, v) and graph.degree(u) == 2 and graph.degree(v) == 2:
            graph.add_edge(u, v)
            r += 1
            if r in r_list:
                improved_graphs[r] = graph.copy()

    return improved_graphs



###### simulate random points ######

def base_ring_from_polygon(
        polygon: Polygon, num_points: int, n_sources: int, rng: np.random.Generator, time_limit: float=None
) -> nx.Graph:
    """
    Generate a base cycle (ring) network from a geographic polygon.
    
    Generates random points within the polygon, solves TSP to form a cycle,
    and randomly selects source nodes.
    
    Parameters
    ----------
    polygon : Polygon
        Geographic region to sample from.
    num_points : int
        Number of nodes in the ring.
    n_sources : int
        Number of source/supply nodes to randomly select.
    time_limit : float, optional
        Time limit (seconds) for TSP solving.
    
    Returns
    -------
    nx.Graph
        Ring graph with 'pos' and 'sources' attributes.
    
    Raises
    ------
    ValueError
        If n_sources < 1.
    """
    if n_sources < 1:
        raise ValueError(f"n_sources should be >= 1")
    random_points = _generate_random_points(polygon, num_points, rng)
    base_ring = _base_ring_from_points(random_points, time_limit=time_limit)
    # Use rng.choice with replacement to mimic `random.choices` behaviour
    sources = list(rng.choice(list(base_ring.nodes), size=n_sources, replace=True))
    base_ring.graph["sources"] = sources
    return base_ring

def _generate_random_points(
        polygon: Polygon,
         num_points: int,
         rng: np.random.Generator,
) -> np.array:
    """
    Generate uniformly random points within a polygon boundary.
    
    Oversamples points in the bounding box and filters to those contained
    in the polygon to achieve efficient rejection sampling.
    
    Parameters
    ----------
    polygon : Polygon
        Shapely Polygon to constrain points to.
    num_points : int
        Target number of points to generate.
    
    Returns
    -------
    np.array
        (num_points, 2) array of (x, y) coordinates.
    """
    # Extract polygon bounds
    min_x, min_y, max_x, max_y = polygon.bounds
    # Generate random points within the bounding box using provided RNG
    factor = 2
    final_points = []
    while len(final_points) < num_points:
        random_x = rng.uniform(min_x, max_x, num_points * factor)  # Oversample
        random_y = rng.uniform(min_y, max_y, num_points * factor)
        random_points = gpd.GeoSeries(gpd.points_from_xy(random_x, random_y))
        contained_points = random_points[random_points.within(polygon)]
        final_points.extend([(p.x, p.y) for p in contained_points])
    return np.array(final_points[:num_points])


####### create base tree and ring ####

def _minimum_spanning_tree(points: np.array, rng: np.random.Generator) -> nx.Graph:
    """Build a minimum spanning tree and attach a deterministic random source.

    Parameters
    ----------
    points : np.array
        (n_points, 2) array of (x, y) coordinates.
    rng : np.random.Generator
        RNG used to pick the source node deterministically.
    """
    dist = squareform(pdist(points))
    G = nx.from_numpy_array(dist)
    T = nx.minimum_spanning_tree(G, weight="weight")
    T.graph["sources"] = [int(rng.choice(list(T.nodes)))]
    nx.set_node_attributes(T, dict(zip(T.nodes, points)), "pos")
    return T




def _base_ring_from_points(points: np.array, time_limit: float=None) -> nx.Graph:
    """
    Create a cycle (ring) graph by solving the Traveling Salesman Problem.
    
    Computes pairwise distances, solves TSP to find an approximate tour,
    converts the tour into a cyclic NetworkX graph, and assigns node positions.
    
    Parameters
    ----------
    points : np.array
        (n_points, 2) array of (x, y) coordinates.
    time_limit : float, optional
        Time limit (seconds) for TSP heuristic solver.
    
    Returns
    -------
    nx.Graph
        Cycle graph with 'pos' attributes set to input coordinates.
    
    """
    # Compute pairwise distances and solve TSP
    dist = squareform(pdist(points))
    cycle, cost = solve_tsp_lin_kernighan(dist, verbose=False, time_limit=time_limit)
    # Convert TSP tour (node sequence) to graph edges
    cycle_edges  = zip(cycle, np.roll(cycle, 1))
    pos = dict(enumerate(points))
    base_ring = nx.from_edgelist(cycle_edges)
    nx.set_node_attributes(base_ring, pos, "pos")
    return base_ring



