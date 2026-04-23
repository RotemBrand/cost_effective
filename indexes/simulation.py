import numpy as np
import networkx as nx
import pandas as pd
from tqdm import tqdm
import heapq
from typing import NamedTuple, Literal

type Event = tuple[float, int, bool]
type RelType = Literal["saidi", "pairwise"]

class RelSimulationResult(NamedTuple):
    rel_result: float
    times: np.ndarray
    cum_rel: np.ndarray
    rel_type: RelType


def _key_edge(e: tuple) -> tuple:
    return tuple(sorted(e))

def _key_edges(e_list: set[tuple]) -> tuple[tuple]:
    return tuple(map(_key_edge, e_list))

class ConnectedComponents:
    def __init__(self, graph: nx.Graph, source: int | None, weight_attr: str, rel_type: RelType):
        # validate
        if rel_type not in ["saidi", "pairwise"]:
            raise ValueError(f"rel_type should be saidi or pairwise not {rel_type}")
        if rel_type == "saidi" and (source is None):
            raise ValueError("source must be not None in rel_type='saidi'")
        if (source is not None) and (source not in graph):
            raise ValueError(f"source {source} not in graph")
        if rel_type == "pairwise":
            print(f"Warning! pairwise rel currently not support weighted nodes")

        # init                     
        self.graph = nx.Graph(graph)
        self.source = source
        self.rel_type: RelType = rel_type
        self.down_edges = set()
        self.lru_cash = {}
        self.weight_attr = weight_attr
        self.total_weight = sum(data[weight_attr] for _, data in graph.nodes.items())        

    def add_edge(self, e: tuple):
        if self.graph.has_edge(*e):
            raise ValueError("Try to add existing edge")
        self.graph.add_edge(*e)
        if _key_edge(e) not in self.down_edges:
            raise ValueError("Try to remove non existing edge form self.down_edges")
        self.down_edges.remove(_key_edge(e))
    
    def remove_edge(self, e: tuple):
        if not self.graph.has_edge(*e):
            raise ValueError("Try to remove non-existing edge")
        self.graph.remove_edge(*e)
        if _key_edge(e) in self.down_edges:
            raise ValueError("Try to add existing edge form self.down_edges")
        self.down_edges.add(_key_edge(e))
    
    def reliability(self) -> float:
        # check if exist in lru cash
        key_down_edges = _key_edges(self.down_edges)
        if key_down_edges in self.lru_cash:
            return self.lru_cash[key_down_edges]

        # calculate and add to cash - saidi
        if self.rel_type == "saidi":
            source_comp = nx.node_connected_component(self.graph, self.source)
            connected_weight = sum(self.graph.nodes[node][self.weight_attr] for node in source_comp)
            disconnected_weight =  1 - connected_weight / self.total_weight
            if len(key_down_edges) <= 2:
                self.lru_cash[key_down_edges] = disconnected_weight
            return disconnected_weight
        
        # calculate and add to cash - pairwise
        elif self.rel_type == "pairwise":
            comps = nx.connected_components(self.graph)
            connected_pairs = 0
            for comp in comps:
                n = len(comp)
                connected_pairs += n * (n - 1) // 2
            n_nodes = len(self.graph)
            disconnected_pairs = 1 - connected_pairs / (n_nodes * (n_nodes - 1) / 2)
            if len(key_down_edges) <= 2:
                self.lru_cash[key_down_edges] = disconnected_pairs
            return disconnected_pairs
        raise ValueError(f"rel_type should be saidi or pairwise")
    

    def reset(self):
        for e in list(self.down_edges):
            self.add_edge(e)



def simulate_rel(
    graph: nx.Graph,
    source,
    rel_type: RelType,
    prob_attr: str,
    weight_attr: str,
    T_days: float,
    mean_cycle_days: float,
    *,
    components: ConnectedComponents | None=None,
    seed: int | None = None,
    rng: np.random.Generator | None = None,
    show_progress: bool = True,
) -> tuple[RelSimulationResult, ConnectedComponents]:
    """
    Event-driven reliability simulation (no dt), NetworkX version, with node weights.

    Returns
    -------
    mean_saidi : float
        Time-average of disconnected *weight fraction* over [0, T_days].
    change_times : np.ndarray
        Times where disconnected weight changed (includes 0.0).
    change_disc_frac : np.ndarray
        Disconnected weight fraction at each change time.
    """
    if rng is None:
        rng = np.random.default_rng(seed)

    if T_days <= 0:
        raise ValueError("T_days must be > 0")
    if mean_cycle_days <= 0:
        raise ValueError("mean_cycle_days must be > 0")

    # ---- node weights ----
    node_weight = {}
    for v in graph.nodes():
        if weight_attr not in graph.nodes[v]:
            raise KeyError(f"Missing node attribute '{weight_attr}' on node {v}")
        node_weight[v] = float(graph.nodes[v][weight_attr])


    edges = list(graph.edges())
    m = len(edges)

    # ---- rates (per day) ----
    alpha = 1.0 / float(mean_cycle_days)

    p = np.empty(m, dtype=float)
    for i, (u, v) in enumerate(edges):
        if prob_attr not in graph[u][v]:
            raise KeyError(f"Missing edge attribute '{prob_attr}' on edge {(u, v)}")
        p_i = float(graph[u][v][prob_attr])
        p[i] = float(np.clip(p_i, 1e-12, 1.0 - 1e-12))

    lam = alpha * p
    rep = alpha * (1.0 - p)

    # ---- live graph ----
    if components is None:
        components = ConnectedComponents(graph, source=source, weight_attr=weight_attr, rel_type=rel_type)
    else:
        if components.source != source:
            raise ValueError(f"Passed components object with source {components.source}, which is different then {source}")
        if components.rel_type != rel_type:
            raise ValueError(f"Passed components object with rel_type {components.rel_type}, which is different then {rel_type}")
        components.reset()



    # ---- initial state ----
    reliability = components.reliability()

    change_times = [0.0]
    change_vals = [reliability]

    # ---- generate all events ----
    events: list[tuple[float, int, bool]] = []

    for i in range(m):
        t = 0.0
        state_up = True
        while True:
            rate = lam[i] if state_up else rep[i]
            t += float(rng.exponential(1.0 / rate))
            if t > T_days:
                break
            state_up = not state_up
            events.append((t, i, state_up))

    events.sort(key=lambda x: x[0])

    # ---- main loop ----
    area = 0.0   # integrates disconnected_weight(t)
    t_prev = 0.0

    iterator = tqdm(events, total=len(events), desc="Events") if show_progress else events
    for t_evt, i, new_up in iterator:
        # integrate
        area += (t_evt - t_prev) * reliability
        t_prev = t_evt

        # apply edge toggle
        if new_up:
            components.add_edge(edges[i])
        else:
            components.remove_edge(edges[i])

        # recompute disconnected weight
        new_reliability = components.reliability()
        if new_reliability != reliability:
            reliability = new_reliability
            change_times.append(t_evt)
            change_vals.append(reliability)

    # integrate remainder
    if t_prev < T_days:
        area += (T_days - t_prev) * reliability
        change_times.append(T_days)
        change_vals.append(reliability)

    mean_saidi = area / (T_days)

    # return the cumulative saidi
    times, cumulative_saidi = compute_cumulative_rel(np.asarray(change_times), np.asarray(change_vals))
    return RelSimulationResult(mean_saidi, times, cumulative_saidi, rel_type), components





def compute_cumulative_rel(change_times, change_disc_frac):
    """
    Compute cumulative rel(t) from piecewise-constant disconnected fraction.

    Parameters
    ----------
    change_times : array-like, shape (K,)
        Times where the disconnected fraction changes.
        Must start at 0 and end at T.
    change_disc_frac : array-like, shape (K,)
        Disconnected fraction on [change_times[i], change_times[i+1]).

    Returns
    -------
    times : np.ndarray, shape (K,)
        Same as change_times.
    cumulative_saidi : np.ndarray, shape (K,)
        Cumulative SAIDI evaluated at each time.
        cumulative_saidi[i] = (1 / times[i]) * ∫₀^{times[i]} d(t) dt
        with cumulative_saidi[0] = 0.
    """
    t = np.asarray(change_times, dtype=float)
    d = np.asarray(change_disc_frac, dtype=float)

    if t.ndim != 1 or d.ndim != 1:
        raise ValueError("Inputs must be 1D arrays")
    if len(t) != len(d):
        raise ValueError("change_times and change_disc_frac must have same length")
    if t[0] != 0:
        raise ValueError("change_times must start at 0")
    if np.any(np.diff(t) < 0):
        raise ValueError("change_times must be non-decreasing")

    # interval lengths
    dt = np.diff(t)

    # cumulative integral of disconnected fraction
    area = np.zeros_like(t)
    area[1:] = np.cumsum(dt * d[:-1])

    # cumulative SAIDI(t) = area(t) / t
    cumulative_rel = np.zeros_like(t)
    mask = t > 0
    cumulative_rel[mask] = area[mask] / t[mask]

    return t, cumulative_rel



def simulate_reliability(
    graph: nx.Graph,
    source,
    mean_time_to_failure: float,
    mean_time_to_repair: float,
    T: float,
    dt: float,
    seed: int | None = None,
    show_progress: bool = True,
    rng: np.random.Generator = None,
):
    """
    Faster SAIDI simulation:
      - Maintains a single mutable graph G and updates it only when events occur.
      - Uses a heap of (event_time, edge_index, event_type) to avoid scanning arrays.
      - Computes connected component only if something changed since last sample time.
    Returns:
      times (np.ndarray), saidi_vector (np.ndarray)
    """
    # RNG and rates
    if rng is None:
        rng = np.random.default_rng(seed)
    lambd = 1.0 / mean_time_to_failure
    repar = 1.0 / mean_time_to_repair

    # Static data
    edges = list(graph.edges())
    n_edges = len(edges)
    n_nodes = graph.number_of_nodes()
    times = np.arange(0.0, T, dt)

    # Working state + live graph
    edge_states = np.ones(n_edges, dtype=bool)  # True=up, False=down
    G = nx.Graph()
    G.add_nodes_from(graph.nodes())
    G.add_edges_from(edges)  # initially all up

    # Build a heap of next events: (time, i, etype) where etype: 0=fail, 1=repair
    # Start by scheduling a failure for each edge
    heap = []
    first_event = rng.exponential(1.0 / lambd, size=n_edges)
    for i in range(n_edges):
        # t_event = rng.exponential(1.0 / lambd)
        t_event = first_event[i]
        heap.append((t_event, i, 0))
    heapq.heapify(heap)

    saidi_vector = []
    progress = tqdm if show_progress else (lambda x, **k: x)
    last_component_size = None   # cache last |connected to source|
    changed_since_last_check = True  # force initial connectivity check

    # Small epsilon to treat float comparison robustly
    EPS = 1e-12
    comp = nx.node_connected_component(G, source)
    for t in progress(times):
        # Process all events up to time t (inclusive with a tiny epsilon)
        step_changed = False
        while heap and heap[0][0] <= t + EPS:
            t_event, i, etype = heapq.heappop(heap)
            u, v = edges[i]

            if etype == 0:  # failure
                if edge_states[i]:
                    edge_states[i] = False
                    # Update the live graph only if edge was present
                    if G.has_edge(u, v):
                        G.remove_edge(u, v)
                    step_changed = True
                # schedule its repair
                heapq.heappush(heap, (t_event + rng.exponential(1.0 / repar), i, 1))

            else:  # repair
                if not edge_states[i]:
                    edge_states[i] = True
                    # Re-add only if not already present
                    if not G.has_edge(u, v):
                        G.add_edge(u, v)
                    step_changed = True
                # schedule next failure
                heapq.heappush(heap, (t_event + rng.exponential(1.0 / lambd), i, 0))

        # Only recompute connectivity if something changed (now or earlier not checked)
        if step_changed or changed_since_last_check or last_component_size is None:
            if source in G:
                # BFS/CC once
                comp = nx.node_connected_component(G, source)
                connected_size = len(comp)
            else:
                connected_size = 0
            last_component_size = connected_size
            changed_since_last_check = False
        else:
            # nothing changed; reuse last_component_size
            connected_size = last_component_size

        n_disconnected = n_nodes - connected_size
        saidi_vector.append(n_disconnected / (n_nodes - 1))

        # If there were events after the last connectivity check but before next t,
        # the flag will be set on next loop. We set it here only if you'd like to
        # force a re-check next time; in our design we already checked when changed.
        # (No action needed.)
    return times, np.array(saidi_vector)


def simulate_reliability_multi_p(
    graph: nx.Graph,
    p_list: list,
    source: int,
    mean_time_to_repair: float=0.5,
    T: float=365 * 5,
    dt: float=0.1,
    seed: int=None,
    file_name: str=None,
    show_progress_sim: bool=True,
    rng: np.random.Generator = None,
):
    # Collect frames and concat once to avoid repeated-concat issues and warnings
    frames: list[pd.DataFrame] = []
    for p in p_list:
        if p > 1e-9:
            mean_time_to_failure = mean_time_to_repair * (1 - p) / p
        else:
            mean_time_to_failure = T
        times, saidi = simulate_reliability(
            graph=graph,
            source=source,
            mean_time_to_failure=mean_time_to_failure,
            mean_time_to_repair=mean_time_to_repair,
            T=T,
            dt=dt,
            seed=seed,
            show_progress=show_progress_sim,
            rng=rng,
        )
        temp_data = pd.DataFrame({
            't': times,
            'saidi': saidi,
            'p': [p for _ in range(len(times))]
        })
        # ensure consistent column order and drop any all-NA columns (defensive)
        temp_data = temp_data[['t', 'p', 'saidi']].copy()
        frames.append(temp_data)
    # Concatenate once (handles empty frames list safely)
    if len(frames) > 0:
        data = pd.concat(frames, ignore_index=True, sort=False)
    else:
        data = pd.DataFrame(columns=['t', 'p', 'saidi'])
    if file_name is not None:
        data.to_csv(file_name)
    return data

def saidi_using_simulation(G: nx.Graph, p_list: list[float], sources: list, seed: int=10, show_progress=False) -> dict[float, float]:
    """
    Estimate SAIDI for multiple failure probabilities using Monte-Carlo simulation.
    
    Contracts all secondary sources into the primary source node, then runs
    the stochastic reliability simulator for each p value.
    
    Parameters
    ----------
    G : nx.Graph
        Network graph to simulate.
    p_list : list[float]
        List of failure probabilities.
    sources : list
        List of source nodes (first is primary).
    seed : int, optional
        Random seed for simulator.
    
    Returns
    -------
    dict[float, float]
        Mapping from p value to estimated mean SAIDI.
    """
    # Contract sources into primary node to unify multiple sources
    G_copy = nx.MultiGraph(G)
    for source in sources[1:]:
        nx.contracted_nodes(G_copy, sources[0], source, self_loops=False, copy=False)
    # make simulation
    print(f"TODO!!! test with multigraph")
    sim_res = simulate_reliability_multi_p(
        G_copy,
        p_list,
        source=sources[0],
        seed=seed,
        mean_time_to_repair=0.5,
        T=365*5,
        dt=0.05,
        show_progress_sim=show_progress
    )
    saidi_by_p = dict(sim_res.groupby('p')['saidi'].mean())
    return saidi_by_p

