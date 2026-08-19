import numpy as np
import networkx as nx
import gurobipy as gp
from gurobipy import GRB
from scipy.spatial import Voronoi, cKDTree
from itertools import combinations
from .utilities import toc, ILPNoSolutionError 
from typing import List

type Node = int
type Chain = int
type Trip = tuple[Node, Node, Chain]


def get_optimal_strc_trips(
        centers: np.array,
        points: np.array,
        n_nodes: int,
        strc_n_init_iters: int,
        exact_vertices: bool,
        debug: bool,
        source_node: int,
        max_trip_vertices_per_center: int | None = None,
    ) -> List[Trip]:
    """Compute a minimal cubic, 3-edge-connected structure expressed as trips.

    This function is the second stage of the `optimal_network_from_points`
    pipeline. The algorithm works in three conceptual phases:

    1. Cluster the input points (`points`) into chains; the cluster centers
       are supplied here as `centers` (result of the k-means stage).
    2. Build a dense "trip" graph where each trip represents a connection
       between two original points (nodes) going through a chain center.
       Concretely, a Trip is a tuple ``(u, v, c)`` meaning the path from
       point u to point v via center c; its weight is the sum of Euclidean
       distances from u to c and from c to v.
    3. Solve an ILP to choose a minimal-weight subgraph of these trips that is
       cubic (degree 3 at selected vertices), 3-edge-connected, and respects
       constraints (one trip per center, fixed number of structure nodes).

    Parameters
    ----------
    centers : np.array
        Array of cluster centers (chain centers) produced by the clustering
        stage; shape (k, 2).
    points : np.array
        Array of original point coordinates; shape (n_points, 2).
    n_nodes : int
        Number of structure (vertex) nodes to select in the resulting
        3-connected structure.
    strc_n_init_iters : int
        Number of ILP initialization iterations / cutting-plane iterations to
        perform before switching strategies (controls runtime vs quality).
    exact_vertices : bool
        If True, the ILP is formulated over exact candidate vertices; if
        False a floating-vertex model is used first and then refined.
    debug : bool
        Enable debug printing and progress information.
    source_node : int
        Index of the original point that should be forced into the structure
        (if not None). This ensures the chosen structure includes the source.

    Returns
    -------
    List[Trip]
        A list of chosen trips; each trip is (point_idx_u, point_idx_v, center_idx)
        where the point indices refer to positions in the original `points`
        array and center_idx identifies which chain center the trip passes through.

    Raises
    ------
    ILPNoSolutionError
        If no feasible cubic 3-edge-connected structure can be found under the
        specified constraints.
    """

    # check for ring graph
    if n_nodes == 1:
        return  [[0, 0, 0]]
    # construct a graph of all the optinal vertices-chains trips
    toc(debug=debug)
    max_n_vertices = n_nodes if exact_vertices else None
    strc = _build_voronoi_trip_graph(
        centers,
        points,
        min_n_vertices=n_nodes,
        max_n_vertices=max_n_vertices,
        source_node=source_node,
        max_trip_vertices_per_center=max_trip_vertices_per_center,
    )
    # use the trips to construct a 3 connected structure that pass all the chains
    toc("build_voronoi_trip_graph", debug)
    chosen_trips = build_cubic_3edge_connected(
        strc,
        n_nodes = n_nodes,
        n_init_iters=strc_n_init_iters,
        exact_vertices=exact_vertices,
        source_node=source_node,
        debug=debug
    )
    toc("build_cubic_3edge_connected", debug)
    return chosen_trips

##### build strc #######

def _build_voronoi_trip_graph(
        centers: np.array,
        points: np.array,
        min_n_vertices: int,
        source_node: int,
        max_n_vertices: int=None,
        seed: int=42,
        max_trip_vertices_per_center: int | None = None,
    ) -> nx.MultiGraph:
    # calculate Voronoi cells
    vor = Voronoi(centers)
    centers_dict = dict(enumerate(centers))

    # filter only inner vertices
    vertices = vor.vertices
    vertices = vertices[vertices[:, 0] <= centers[:, 0].max()]
    vertices = vertices[vertices[:, 0] >= centers[:, 0].min()]
    vertices = vertices[vertices[:, 1] <= centers[:, 1].max()]
    vertices = vertices[vertices[:, 1] >= centers[:, 1].min()]

    # add midpoints if there are not enougth vertices
    if len(vertices) < min_n_vertices:
        rp = vor.ridge_points
        midpoints = (points[rp[:, 0]] + points[rp[:, 1]]) / 2
        # choose random midpoints
        rng = np.random.default_rng(0)
        idx = rng.choice(list(range(len(midpoints))), min_n_vertices - len(vertices))
        midpoints = midpoints[idx]
        vertices = np.r_[vertices, midpoints]

    # choose random vertices if there are too much
    if max_n_vertices and len(vertices) > max_n_vertices:
        rng = np.random.default_rng(seed)
        vertices_idx = rng.choice(list(range(len(vertices))), max_n_vertices)
        vertices = vertices[vertices_idx]

    # find closest point for each vertex
    chosen_points_id = _closest_point_to_each_vertex(points=points, vertices=vertices)
    if (source_node is not None) and (source_node not in chosen_points_id):
        chosen_points_id = np.concatenate([chosen_points_id, [source_node]])
        if max_n_vertices and (len(chosen_points_id) > max_n_vertices):
            rng = np.random.default_rng(seed)
            id_to_remove = set(rng.choice(chosen_points_id, len(chosen_points_id) - max_n_vertices))
            chosen_points_id = np.array([x for x in chosen_points_id if x not in id_to_remove])
    chosen_points = points[chosen_points_id]


    trips_graph = nx.MultiGraph()
    for i, pos in enumerate(chosen_points):
        trips_graph.add_node(i, pos=pos, point_idx=chosen_points_id[i])

    if max_trip_vertices_per_center is not None:
        if max_trip_vertices_per_center < 2:
            raise ValueError("max_trip_vertices_per_center must be at least 2")
        k = min(max_trip_vertices_per_center, len(chosen_points))
        tree = cKDTree(chosen_points)
        for ck, c in centers_dict.items():
            _, nearest = tree.query(c, k=k)
            nearest = np.atleast_1d(nearest).astype(int).tolist()
            for i, j in combinations(nearest, 2):
                pos_i, pos_j = chosen_points[i], chosen_points[j]
                weight = np.linalg.norm(pos_i - c) + np.linalg.norm(pos_j - c)
                trips_graph.add_edge(i, j, ck, cpos=c, weight=weight)
        return trips_graph

    # make a graph of all v1 -> c -> v2 trips
    for (i, j) in combinations(list(range(len(chosen_points_id))), 2):
        pos_i, pos_j = chosen_points[i], chosen_points[j]
        for ck, c in centers_dict.items():
            weight = np.linalg.norm(pos_i - c) + np.linalg.norm(pos_j - c)
            trips_graph.add_edge(i, j, ck, cpos=c, weight=weight)
    return trips_graph
    

def _closest_point_to_each_vertex(points: np.array, vertices: np.array) -> np.array:
    """Return an array of points_id of the closest point per vertex with no duplicates"""
    tree = cKDTree(points)

    points_idx = [] #point_id
    selected_points_idx = set()
    for vertex_id, vertex in enumerate(vertices):
        k = 1
        while (k == 1) or (point_id in selected_points_idx): 
            _, (point_id, ) = tree.query(vertex, k=[k])
            k += 1
        selected_points_idx.add(point_id)
        points_idx.append(point_id)
    return np.array(points_idx)


########## find 3 connected ##########

def key_edges(edges):
    return tuple(sorted([key(e) for e in edges]))

def key(e):
    if len(e) == 3:
        return (*sorted(e[:2]), e[2])
    return tuple(sorted(e))

def _find_trip_source_node(trips_graph: nx.MultiGraph, source_node: int) -> int:
    if source_node is None:
        return None
    for node, point_idx in nx.get_node_attributes(trips_graph, "point_idx").items():
        if point_idx == source_node:
            return node
    return None

def build_cubic_3edge_connected(
        trips_graph: nx.MultiGraph,
        n_nodes: int,
        n_init_iters: int=2,
        exact_vertices: bool=False,
        source_node: int=None,
        debug: bool=False
    ) -> List[Trip]:
    trip_source = _find_trip_source_node(trips_graph, source_node)
    if exact_vertices:
        m, x, edges = _init_exact_vertices_model(trips_graph, n_nodes)
    else: 
        m, x, y, edges, nodes = _init_floating_vertices_model(trips_graph, n_nodes, source_node=trip_source)
    # Iterative solve
    it = 0
    while True:
        m.optimize()
        if m.Status != GRB.OPTIMAL: break
        status = _add_edge_connectivity_cuts(m, x, edges, n_nodes, debug=debug)
        chosen_edges = [e for e in edges if x[e].X > 0.5]

        if not status:
            break
        if it >= n_init_iters:
            # igonore existing edges in teh objective to speed the converge
            chosen_edges = [e for e in edges if x[e].X > 0.5]
            m.setObjective(gp.quicksum(trips_graph.edges[e]['weight'] * x[e] for e in edges if e not in chosen_edges), GRB.MINIMIZE)
            
            # transform a floting model into exact vetrices model
            if (not exact_vertices) and it == n_init_iters:
                chosen_nodes = [v for v in nodes if y[v].X > 0.5]
                subgraph = nx.subgraph(trips_graph, chosen_nodes).copy()
                return build_cubic_3edge_connected(
                    trips_graph=subgraph,
                    n_nodes=n_nodes,
                    n_init_iters=1,
                    exact_vertices=True,
                    debug=debug
                )
        it += 1

    # handle non fessible
    status = m.Status
    if status != GRB.Status.OPTIMAL:
        raise ILPNoSolutionError("Can not find 3 regular structure")
    
    # extract trips from solution
    chosen_edges = [e for e in edges if x[e].X > 0.5]
    H = nx.Graph()
    H.add_edges_from([e[:2] for e in chosen_edges])
    chosen_trips: List[Trip] = [
        (trips_graph.nodes[u]["point_idx"], trips_graph.nodes[v]["point_idx"], c)
        for u, v, c in chosen_edges
    ]
    return chosen_trips


def _init_floating_vertices_model(
    trips_graph: nx.MultiGraph,
    n_nodes: int,
    source_node: int=None
):
    # init model
    m = gp.Model("Cubic3EdgeConn")
    m.Params.OutputFlag = 0

    # Decision vars
    edges = list(map(key, trips_graph.edges))
    nodes = list(trips_graph.nodes)
    x = m.addVars(edges, vtype=GRB.BINARY, name="x")  # edges
    y = m.addVars(nodes, vtype=GRB.BINARY, name="y")  # nodes

    # 1. Select exactly n_nodes and n_edges
    n_edges = int(1.5 * n_nodes)
    m.addConstr(gp.quicksum(y[v] for v in nodes) == n_nodes, name="node_count")
    m.addConstr(gp.quicksum(x[e] for e in edges) == n_edges, name="edge_count")

    # 2. Degree 3 constraint
    for v in nodes:
        m.addConstr(gp.quicksum(x[key(e)] for e in trips_graph.edges(v, keys=True)) == 3*y[v], name=f"deg_{v}")

    # 3. Each center appears exactly once
    centers = {}
    for u, v, c in edges:
        centers.setdefault(c, []).append((u,v, c))
    for c, c_edges in centers.items():
        m.addConstr(gp.quicksum(x[e] for e in c_edges) == 1, name=f"center_{c}")

    # 4. Edge-node consistency
    for e in edges:
        m.addConstr(x[e] <= y[e[0]])
        m.addConstr(x[e] <= y[e[1]])

    # 5. choose at most one key from each multiedge
    if n_nodes >= 3:
        for u, v in nx.Graph(trips_graph).edges:
            m.addConstr(gp.quicksum(x[key((u, v, c))] for c in trips_graph[u][v]) <= 1)

    # 6. force to choose source node
    if source_node is not None:
        m.addConstr(y[source_node] == 1)

    # Objective: minimize sum of weights
    m.setObjective(gp.quicksum(trips_graph.edges[e]['weight'] * x[e] for e in edges), GRB.MINIMIZE)

    return m, x, y, edges, nodes

def _init_exact_vertices_model(
    trips_graph: nx.MultiGraph,
    n_nodes: int,
):
    # init model
    m = gp.Model("Cubic3EdgeConn")
    m.Params.OutputFlag = 0

    # Decision vars
    edges = list(map(key, trips_graph.edges))
    nodes = list(trips_graph.nodes)
    x = m.addVars(edges, vtype=GRB.BINARY, name="x")  # edges

    # 1. Select exactly n_edges
    n_edges = int(1.5 * n_nodes)
    m.addConstr(gp.quicksum(x[e] for e in edges) == n_edges, name="edge_count")

    # 2. Degree 3 constraint
    for v in nodes:
        m.addConstr(gp.quicksum(x[key(e)] for e in trips_graph.edges(v, keys=True)) == 3, name=f"deg_{v}")

    # 3. Each center appears exactly once
    centers = {}
    for u, v, c in edges:
        centers.setdefault(c, []).append((u, v, c))
    for c, c_edges in centers.items():
        m.addConstr(gp.quicksum(x[e] for e in c_edges) == 1, name=f"center_{c}")

    # 5. choose at most one key from each multiedge
    if n_nodes >= 3:
        for u, v in nx.Graph(trips_graph).edges:
            m.addConstr(gp.quicksum(x[key((u, v, c))] for c in trips_graph[u][v]) <= 1)


    # Objective: minimize sum of weights
    m.setObjective(gp.quicksum(trips_graph.edges[e]['weight'] * x[e] for e in edges), GRB.MINIMIZE)

    return m, x, edges

# Cutting-plane loop for 3-edge-connectivity
def _add_edge_connectivity_cuts(m, x, edges, n_nodes, debug: bool):
    # Extract current solution
    H = nx.Graph()
    chosen_edges = [e for e in edges if x[e].X > 0.5]
    H.add_edges_from([e[:2] for e in chosen_edges])

    # r = 2
    if n_nodes == 2:
        return False

    if len(H) == 0: return False
    comps = list(nx.k_edge_components(H, 3))
    if debug:
        print(f"len(comps) = {len(comps)}")
    if len(comps) == 1:
        return False  # no violated cuts

    # Find specific minimum cut
    for comp in comps:
        # Add constraint: sum of these edges >= 3
        m.addConstr(gp.quicksum(x[e] for e in edges if len(set(e[:2]) & comp) == 1) >= 3)
    return True

