"""
The main file of the GraphRel object used to calculate the graph reliability
"""

import networkx as nx
import numpy as np
from enum import Enum
from copy import copy
from typing import Dict, List, Tuple, Union, Literal
from numpy.typing import ArrayLike
from functools import reduce
from itertools import combinations
from abc import ABCMeta
import indexes.tree_decomposition as TD
import indexes.utilities as utilities
import indexes.probs as PROBS
from indexes.simulation import simulate_rel, RelSimulationResult
from tqdm import tqdm


multiply = lambda x, y: x * y
add = lambda x, y: x + y
_SUPPORTED_REL_KINDS_WITH_STRUCTURE = ["saidi"]
Method = Literal["td", "pivot"]
type RelType = Literal["saidi", "pairwise"]

def array_to_size(arr: np.array, s: int) -> np.array:
    """## Resize arr to length s. if len(arr) < s, return only
    the first s elements. Otherwise, pad with zeros

    ### Args:
        - `arr (np.array)`: array
        - `s (int)`: new size for the array

    ### Returns:
        - `np.array`: New array from length s
    """
    if not hasattr(arr, "__len__"):
        input_arr = np.array([arr])
    else:
        input_arr = np.array(arr)
    if s <= len(input_arr):
        return input_arr[:s]
    res = np.zeros(s)
    res[: len(input_arr)] = input_arr
    return res


def get_degree(poly_list: List[PROBS.Poly]) -> int:
    """Get a list of polynoms and output their polynom degree"""
    poly_list_trim = [np.trim_zeros(poly.prob, "b") for poly in poly_list]
    return sum([len(poly) - 1 for poly in poly_list_trim])


class RelKind(Enum):
    """An enum class that stores the reliability indexes types
    - SAIDI: the mean number of disconnected nodes from the source
    - ALL_TERMINAL: The probability that the graph is disconnected
    - PAIRWISE: The mean number of disconnected nodes pairs
    - EDGES_REL: The reaching probabilities to each edge
    - NODES_REL: THe disconnected probability of each node
    """

    SAIDI = "saidi"
    ALL_TERMINAL = "all_terminal"
    PAIRWISE = "pairwise"
    EDGES_REL = "edges_rel"
    NODES_REL = "nodes_rel"

    @classmethod
    def valid_values(cls):
        """Return all the valid class values"""
        return list(cls.__members__.keys())


class GraphRel:
    r"""
    The class contain all the information to calculate the
    graph reliability. Each object contain a graph and related configuration used to
    calculate its reliability.

    __Attributes__:
        * original_graph(nx.MultyGraph): The original graph used to define the object
        * graph(nx.MultyGraph): The basis graph used to calculate. Contains one source after contraction,
            and the nodes and edges attributes. Also, if original graph contain non-integers nodes, its relable it.
        * prob_class (PROBS.Prob): The probability type of the edges failure probability (Float/Array/Poly...).
        * source(any): The source of self.graph after the sources contraction.
        * max_fail (int): The maximal number of failing edges allowed for the computation. Used to approximate the reliability for
            small failure probability. Can be not None only if prob_type in Probs.Poly
        * recursive_count (int): The number of times that culc_rel_recursive was activated.
        * res (Dict): The result of the reliability calculations
    """

    def __init__(
        self,
        graph: nx.MultiGraph | nx.Graph,
        nodes_weight: Union[Dict, float] | None = None,
        edges_prob: Union[Dict, PROBS.Prob] | None = None,
        edges_weight: Union[Dict, float] | None = None,
        max_fail: int| None = None,
        sources: ArrayLike = [0],
    ):
        """
        Initialize the GraphRel object for reliability calculations.

        This method sets up the original graph and creates a transformed graph
        with attributes such as node weights and edge failure probabilities.
        If there are multiple sources, it contracts them into a single source.
        Non-integer nodes are relabeled for faster computation.

        Args:
            graph (nx.MultiGraph): The input graph used for reliability calculations.
            nodes_weight (Union[Dict, float], optional): Node weights as a dictionary 
                or scalar. Defaults to None, which sets all weights to 1.
            edges_prob (Union[Dict, PROBS.Prob], optional): Edge failure probabilities 
                as a dictionary or scalar. Defaults to None, which sets probabilities 
                to default values.
            edges_weight (Union[Dict, float], optional): Edge weights as a dictionary
                or scalar. Defaults to None, which sets all edge weights to 0.
            max_fail (int, optional): Maximum number of failing edges allowed. Only 
                applicable when the probability type supports it. Defaults to None.
            sources (ArrayLike, optional): List of source nodes for the graph. 
                Defaults to [0].
        """
        # init self.graph
        self.max_fail = max_fail
        self.original_graph = graph.copy()
        self.prob_class = None
        self.graph = self.init_graph_attributes(
            copy(graph), copy(nodes_weight), copy(edges_weight), copy(edges_prob)
        )
        # set the source nodes
        for node in sources:
            if node not in self.original_graph.nodes:
                raise ValueError(f"Source {node} not in graph")
        sources_int = [self.nodes_to_int[source] for source in sources]
        self.graph = contracted_nodes_list(self.graph, sources_int)
        self.source = sources_int[0]
        # get only the connected component of the source
        self.graph.nodes[self.source]["weight"] = 0
        graph = nx.subgraph(
            self.graph, nx.node_connected_component(self.graph, self.source)
        )
        self.recursive_count = 0
        self._td = None
        self._res = {}

    @property
    def td(self):
        """The tree decomposition of self.graph.

        Returns
        -------
        nx.Graph
            A graph where each node is a graph subset
        """
        if self._td is None:
            self._td = TD.get_td_with_attributes(self)
        return self._td

    def init_graph_attributes(
        self,
        graph: nx.MultiGraph,
        nodes_weight: Union[Dict, float],
        edges_weight: Union[Dict, float],
        edges_prob: Union[Dict, PROBS.Prob],
    ) -> nx.MultiGraph:
        """
        Initialize graph attributes for reliability calculations.

        Adds attributes such as node weights, edge failure probabilities, 
        and other metadata to the graph.

        Args:
            graph (nx.MultiGraph): The input graph to configure.
            nodes_weight (Union[Dict, float]): Node weights as a dictionary or scalar.
            edges_weight (Union[Dict, float]): Edge weights as a dictionary or scalar.
            edges_prob (Union[Dict, PROBS.Prob]): Edge failure probabilities.

        Returns:
            nx.MultiGraph: The graph with configured attributes.
        """
        if graph.number_of_edges() == 0:
            raise ValueError("Graph without edges")
        # set original edges attributes
        if isinstance(graph, nx.MultiGraph):
            original_edges = {e: e for e in graph.edges}
        elif isinstance(graph, nx.Graph):
            original_edges = {e + (0,): e for e in graph.edges}
        else:
            raise ValueError("graph should be a graph or a multigraph")
        graph = nx.MultiGraph(graph)
        nx.set_node_attributes(graph, {v: {v} for v in graph.nodes}, "component")
        nx.set_edge_attributes(graph, original_edges, "original_edge")
        # set nodes weight
        if nodes_weight is None:
            nodes_weight = {v: 1 for v in graph.nodes}
        elif not isinstance(nodes_weight, dict):
            nodes_weight = {v: nodes_weight for v in graph.nodes}
        nx.set_node_attributes(graph, nodes_weight, "weight")
        for node, node_data in graph.nodes.items():
            if "weight" not in node_data:
                raise ValueError(f"node {node} has no weight")
        # set edge weights and tie metadata
        edges_weight = self._init_edge_attr_from_input(
            graph,
            edges_weight,
            attr="edge_weight",
            default=0.0,
        )
        nx.set_edge_attributes(graph, edges_weight, "edge_weight")
        is_tie = self._init_edge_attr_from_input(
            graph,
            None,
            attr="is_tie",
            default=False,
        )
        nx.set_edge_attributes(graph, is_tie, "is_tie")
        # set edges_prob
        if edges_prob is None:
            edges_prob = {edge: PROBS.Poly([0, 1]) for edge in graph.edges}
        elif not isinstance(edges_prob, dict):
            edges_prob = {edge: edges_prob for edge in graph.edges}
        prob_class, edges_prob_new = define_prob_class_and_init_edges_prob(
            edges_prob, self.max_fail
        )
        self.prob_class = prob_class
        edges_prob_has_key_of_edge = len(next(iter(edges_prob_new.keys()))) == 3
        if not edges_prob_has_key_of_edge:
            edges_prob_new = {(k[0], k[1], 0): v for k, v in edges_prob_new.items()}
        nx.set_edge_attributes(graph, edges_prob_new, "prob")
        for edge, edge_data in graph.edges.items():
            if "prob" not in edge_data:
                raise ValueError(f"prob is not set for {edge}")
        if (self.max_fail is not None) and (
            not self.prob_class.need_to_check_zero_events
        ):
            raise ValueError(
                f"There is not option to approximate the reliability using prob_class = {self.prob_class.__name__}. Please set max_fail to None"
            )
        # convet nodes to integers
        if any([not isinstance(node, int) for node in graph.nodes]):
            nodes_to_int = {node: i for i, node in enumerate(graph.nodes)}
        else:
            nodes_to_int = {node: node for node in graph.nodes}
        graph = nx.relabel_nodes(graph, nodes_to_int)
        nx.set_node_attributes(
            graph, {v: k for k, v in nodes_to_int.items()}, "original_node"
        )
        nx.set_node_attributes(self.original_graph, nodes_to_int, "node_int_label")
        self.nodes_to_int = nodes_to_int
        return graph

    @staticmethod
    def _init_edge_attr_from_input(
        graph: nx.MultiGraph,
        values: Union[Dict, float] | None,
        attr: str,
        default,
    ) -> Dict[tuple, object]:
        if values is None:
            return {
                edge: data.get(attr, default)
                for edge, data in graph.edges.items()
            }
        if not isinstance(values, dict):
            return {edge: values for edge in graph.edges}
        if len(values) == 0:
            return {edge: default for edge in graph.edges}
        values_has_key_of_edge = len(next(iter(values.keys()))) == 3
        if values_has_key_of_edge:
            return {
                edge: values.get(edge, graph.edges[edge].get(attr, default))
                for edge in graph.edges
            }
        return {
            edge: values.get(edge[:2], graph.edges[edge].get(attr, default))
            for edge in graph.edges
        }

    def calc_rel(self, rel_kinds: List[str] = None, method: Method="td") -> dict:
        """Calculate the graph reliabilities indexes according to rel_kinds_list

        Parameters
        ----------
        rel_kinds : List[str],
            A list of different reliability methods. All the methods are enums of RelKind Class
            The reliability methods are:
                - saidi: the average number of disconnected nodes
                - all_terminal: the disconnection probability of the graph(NOT IMPLEMENTED),
                - nodes_rel: a dict of the disconnection probability of each node from the source
                - edges_rel: a dict of the reaching probabilities to the edge nodes
                - pairwise: the mean number of disconnected pair of nodes(NOT IMPLEMENTED)
        method : str, optional
            The method used to calculate the reliability:
                - td: use tree-decomposition methods(good for sparse graphs)
                - pivot: the classical pivot method. A very slow method
            by default 'td'

        Returns
        -------
        dict
            each key of the result is a rel_kind and its result

        Raises
        ------
        ValueError
            if the method is not td or pivot
        """
        self._raise_if_edge_weights_for_exact_methods()
        if method == "td":
            return self.culc_rel_td(rel_kinds=rel_kinds)
        if method == "pivot":
            res = {}  # TODO: make it suuport multiple kinds
            for rel_kind in rel_kinds:
                self.kind = RelKind(rel_kind)
                res[rel_kind] = self.culc_rel_pivot()
            return res
        else:
            raise ValueError("methods can be td/pivot")

    @classmethod
    def reliability_array(cls, graph: nx.Graph, sources: list, p: np.array, kind="saidi"):
        if kind != "saidi":
            raise ValueError("Currently GraphRel.reliability_array supported only saidi reliability")
        gr = GraphRel(graph, sources=sources, edges_prob=PROBS.Array(p))
        return np.array(gr.calc_rel([kind])[kind].prob, dtype='float')

    @classmethod
    def reliability_polynomial(cls, graph: nx.Graph, sources: list, kind="saidi", max_fail=None):
        gr = GraphRel(graph, sources=sources, edges_prob=PROBS.Poly([0, 1]), max_fail=max_fail)
        return gr.calc_rel([kind])[kind]

    def culc_saidi_using_skeleton(self, method="td") -> dict:
        """
        Calculate the SAIDI reliability index using the skeleton graph.
        This function is fast for very sparse graphs

        Args:
            method (str, optional): Method for calculation. Currently supports 'td'. 
                Defaults to 'td'.

        Returns:
            dict: A dictionary with the SAIDI reliability index.
        """
        self._raise_if_edge_weights_for_exact_methods()
        # create the skeleton
        gr_skeleton = utilities.create_gr_of_the_skeleton_graph(self)
        # calculate the reliability of the skeleton
        if method == "pivot":
            raise ValueError("pivot method with skeleton is currently unsopoorted")
        elif method == "td":
            rel_skeleton = gr_skeleton.culc_rel_td(rel_kinds=["saidi", "edges_rel"])
        else:
            raise ValueError("valid methods are td or pivot")
        # set edges weight
        for edge, edge_data in gr_skeleton.graph.edges.items():
            edge_weight = sum(
                [
                    self.graph.nodes[node]["weight"]
                    for node in edge_data["nodes"]
                ]
            )
            gr_skeleton.graph.edges[edge]["weight"] = edge_weight
        # calculate the saidi of the original graph
        saidi = self.prob_class([0.0])
        edges_rel = rel_skeleton["edges_rel"]
        # TODO: sum differeny polynoms len
        for edge, edge_data in gr_skeleton.graph.edges.items():
            saidi += (
                (edges_rel[edge][frozenset()] * 1)
                + (
                    edges_rel[edge][frozenset(edge[0:1])]
                    * utilities.path_graph_saidi(edge_data["subgraph"], source=edge[0])
                )
                + (
                    edges_rel[edge][frozenset(edge[1:2])]
                    * utilities.path_graph_saidi(edge_data["subgraph"], source=edge[1])
                )
                + (
                    edges_rel[edge][frozenset(edge[0:2])]
                    * utilities.path_graph_saidi_two_sources(
                        edge_data["subgraph"], source=edge[0]
                    )
                )
            ) * edge_data["weight"]
        saidi += rel_skeleton["saidi"] * gr_skeleton.total_weight()
        total_weight = self.total_weight()
        if total_weight != 0:
            saidi /= total_weight
        else:
            saidi = self.prob_class([0.0])
        self._res["saidi"] = saidi
        return {"saidi": saidi}

    def culc_rel_td(self, rel_kinds: List[RelKind | str]) -> dict:
        """Use tree decomposition method to calculate the graph reliability.
        Its work by decompose the graph to smaller subgraphs, calculate the reliability on each of them and combine the results

        Parameters
        ----------
        rel_kinds : List[RelKind  |  str]
            A list of the reliability indexes

        Returns
        -------
        dict
            A dict where each key is a rel_kind and its results reliability
        """
        self._raise_if_edge_weights_for_exact_methods()
        rel_kinds_ = [RelKind(rk) for rk in rel_kinds]
        # prepare the td for saidi computation
        if (
            (RelKind.SAIDI in rel_kinds_)
            or (RelKind.EDGES_REL in rel_kinds_)
            or (RelKind.NODES_REL in rel_kinds_)
        ):
            # run_out_pc = RelKind.ALL_TERMINAL not in rel_kinds_
            TD.set_the_final_partitions_collection_for_each_bug(self, True)
            # calculate saidi index
            if RelKind.SAIDI in rel_kinds_:
                nodes_probs = TD.get_nodes_disconnection_prob_from_td(self)
                nodes_weights = nx.get_node_attributes(self.graph, "weight")
                if self.total_weight() != 0:
                    saidi = saidi_from_nodes_prob_and_weight(nodes_probs, nodes_weights)
                else:  # handle zero weight
                    edge_example = next(iter(self.graph.edges))
                    saidi = self.graph.edges[edge_example]["prob"] * 0.0
                self._res[RelKind.SAIDI.value] = saidi
            # calculate edge rel
            if RelKind.EDGES_REL in rel_kinds_:
                edges_probs = TD.set_edges_disconnection_prob(self)
                self._res[RelKind.EDGES_REL.value] = edges_rel_in_original_graph(
                    self, edges_probs
                )
                nx.set_edge_attributes(
                    self.original_graph,
                    self._res[RelKind.EDGES_REL.value],
                    RelKind.EDGES_REL.value,
                )
            # calculate nodes_rel
            if RelKind.NODES_REL in rel_kinds_:
                nodes_probs = TD.get_nodes_disconnection_prob_from_td(self)
                nodes_prob_original_graph = rename_nodes_prob_to_original_labels(
                    self, nodes_probs
                )
                self._res[RelKind.NODES_REL.value] = nodes_prob_original_graph
        if RelKind.ALL_TERMINAL in rel_kinds_:
            raise ValueError("Not supported")
        return {kind.value: self._res[kind.value] for kind in rel_kinds_}

    def calc_rel_simulation(
        self ,
        T_days: float=365*5,
        mean_cycle_days: float=0.5,
        rel_type: RelType="saidi",
        *,
        seed: int | None = None,
        rng: np.random.Generator | None = None,
        show_progress: bool = False,
    ) -> RelSimulationResult | tuple[RelSimulationResult]:
        edge_probs = nx.get_edge_attributes(self.graph, "prob")
        if self.prob_class == PROBS.Float:
            edge_probs = {e: [float(p)] for e, p in edge_probs.items()}
        elif self.prob_class == PROBS.Array:
            edge_probs = {e: p.prob for e, p in edge_probs.items()}
        else:
            raise ValueError(f"rel kind when calculate with simulation should only be Float or Array not {self.prob_class}")

        # edge_probs = {e: [p.prob for p in p_arr] for e, p_arr in edge_probs.items()}
        # iterator
        n_of_p = len(next(iter(edge_probs.values())))
        saidi_arr = []
        if show_progress and (self.prob_class==PROBS.Array):
            iterator = tqdm(range(n_of_p), total=n_of_p, desc="Simulate SAIDI")
        else:
            iterator = range(n_of_p)

        components = None
        for i in iterator:
            current_edge_probs = {e: p_arr[i] for e, p_arr in edge_probs.items()}
            graph_copy = self.graph.copy()
            nx.set_edge_attributes(graph_copy, current_edge_probs, "prob")
            saidi_res, components = simulate_rel(
                nx.MultiGraph(graph_copy),
                source=self.source,
                rel_type=rel_type,
                prob_attr="prob",
                weight_attr="weight",
                edge_weight_attr="edge_weight",
                T_days=T_days,
                mean_cycle_days=mean_cycle_days,
                seed=seed,
                rng=rng,
                components=components,
                show_progress=show_progress and (self.prob_class==PROBS.Float)
            )
            saidi_arr.append(saidi_res)
        if self.prob_class == PROBS.Float:
            return saidi_arr[0]
        else:
            return tuple(saidi_arr)


    def culc_rel_pivot(self):  # TODO: make kinds in the input
        """## The main funtion to calculate the graph reliabiltiy. NOT WORK!

        ### Returns:
            - `any`: The reliability of the graph
        """
        self._raise_if_edge_weights_for_exact_methods()
        rel = self.culc_rel_recursive(self.graph.copy(), 0)
        rel = self.prob_class(rel, level=self.max_fail, reverse=True)
        if self.kind == RelKind.SAIDI:
            total_weight = sum(
                [
                    data["weight"]
                    for node, data in self.graph.nodes(data=True)
                    if node != self.source
                ]
            )
            rel = rel / total_weight
        elif self.kind == RelKind.ALL_TERMINAL:
            pass
        print(f"recursive count = {self.recursive_count}")
        return rel

    def culc_rel_recursive(self, graph: nx.MultiGraph, n_failing_edges: int):
        """## The recursive function used to calculate the reliability. The reliability of
        the graph is calcualte recursively using the reliability of the contracted graph
        and the deleted graph.

        ### Args:
            - `graph (nx.MultiGraph)`: The graph to caluclate its reliability
            - `n_failing_edges (int)`: The number of alredy failed edges. Used for stoping condition
            for case that n_failing_edges > self.level

        ### Returns:
            - `any`: The graph reliability
        """
        self.recursive_count += 1
        if n_failing_edges > self.max_fail:
            return 0
        edge = non_self_edge(graph, self.source)
        if edge is None:
            return self.stoping_condition(graph)
        fail_prob = reduce(
            multiply, [data["prob"] for data in graph[edge[0]][edge[1]].values()]
        )
        deleted_graph_rel = self.rel_deleted_graph(graph, edge, n_failing_edges)
        contracted_graph_rel = self.rel_contracted_graph(graph, edge, n_failing_edges)
        return (
            fail_prob * (deleted_graph_rel - contracted_graph_rel)
            + contracted_graph_rel
        )

    def stoping_condition(self, graph: nx.MultiGraph):
        """Determine the stopping condition for the reliability calculation.

        ###Args:
            graph (nx.MultiGraph): The graph being analyzed.

        ###Raises:
            ValueError: If the RelKind is not valid.

        ###Returns:
            float: The calculated stopping condition value.
        """
        match self.kind:
            case RelKind.SAIDI:
                disconnected_weight = sum(
                    [
                        data["weight"]
                        for node, data in graph.nodes(data=True)
                        if node != self.source
                    ]
                )
                return disconnected_weight
            case RelKind.ALL_TERMINAL:
                is_disconnected = float(graph.number_of_nodes() != 1)
                return is_disconnected
            case _:
                raise ValueError(
                    f"self.kind is {self.kind} which is not valid RelKind {list(RelKind.__members__)}"
                )

    def rel_deleted_graph(
        self, graph: nx.MultiGraph, edge: tuple, n_failing_edges: int
    ):
        """Calculate the reliability of the graph with the specified edge deleted.

        ###Args:
            graph (nx.MultiGraph): The graph being analyzed.
            edge (tuple): The edge to delete.
            n_failing_edges (int): The number of already failed edges.

        ###Returns:
            any: The reliability of the graph with the edge deleted.
        """
        u, v = edge[:2]
        edges_to_remove = [(u, v, key, graph[u][v][key]) for key in graph[u][v]]
        graph.remove_edges_from(edges_to_remove)
        deleted_graph_rel = self.culc_rel_recursive(
            graph, n_failing_edges + len(edges_to_remove)
        )
        graph.add_edges_from(edges_to_remove)
        return deleted_graph_rel

    def rel_contracted_graph(
        self, graph: nx.MultiGraph, edge: tuple, n_failing_edges: int
    ):
        """Calculate the reliability of the graph with the specified edge contracted.

        ###Args:
            graph (nx.MultiGraph): The graph being analyzed.
            edge (tuple): The edge to contract.
            n_failing_edges (int): The number of already failed edges.

        ###Returns:
            any: The reliability of the graph with the edge contracted.
        """

        u, v = edge[:2]
        # contract graph
        edges_to_remove = copy(list(graph.edges(v, keys=True, data=True)))
        edges_to_add = [
            (u, node2, max_key(graph, u, node2) + key + 1, data)
            for node1, node2, key, data in edges_to_remove
            if u != node2
        ]
        v_data = graph.nodes[v]
        graph.remove_edges_from(edges_to_remove)
        graph.add_edges_from(edges_to_add)
        graph.remove_node(v)
        # culc reliabiltiy
        contracted_graph_rel = self.culc_rel_recursive(graph, n_failing_edges)
        # restore the original graph
        graph.add_node(v, **v_data)
        graph.remove_edges_from(edges_to_add)
        graph.add_edges_from(edges_to_remove)
        return contracted_graph_rel

    def culc_rel_approximation(self):
        if self.kind != RelKind.SAIDI:
            raise Exception("Implemented only for SAIDI")
        graph = self.multy_graph_to_graph()
        graph_edges = [edge for edge in graph.edges(data=True)]
        rel = 0
        for k in range(1, self.max_fail + 1):
            for edges_to_remove in combinations(graph_edges, k):
                if (
                    sum([data["multiplicity"] for u, v, data in edges_to_remove])
                    > self.max_fail
                ):
                    continue
                edges_to_remove_no_multiplicty = {
                    frozenset(edge[:2]) for edge in edges_to_remove
                }
                graph = self.graph.copy()
                graph.remove_edges_from([e[:2] for e in edges_to_remove])
                source_comp = nx.node_connected_component(graph, self.source)
                non_source_comp = set(graph.nodes) - set(source_comp)
                graph.add_edges_from(edges_to_remove)

                if len(non_source_comp) == 0:
                    continue
                failing_weight = np.sum(
                    [self.graph.nodes[node]["weight"] for node in non_source_comp]
                )
                failing_probability = reduce(
                    multiply, [data["prob"] for (u, v, data) in edges_to_remove]
                )
                sucsess_probability = reduce(
                    multiply,
                    [
                        1 - data["prob"]
                        for (u, v, data) in graph_edges
                        if frozenset([u, v]) not in edges_to_remove_no_multiplicty
                    ],
                )
                sucsess_probability = self.prob_class(
                    sucsess_probability, level=self.poly_level
                )
                rel += failing_probability * failing_weight * sucsess_probability
        return self.prob_class(rel, reverse=True, level=self.poly_level)

    def multy_graph_to_graph(self):
        multy_graph = self.graph.copy()
        graph = nx.Graph(multy_graph)
        graph_edges = list(graph.edges)
        for edge in graph_edges:
            edges_for_same_nodes = multy_graph[edge[0]][edge[1]]
            graph.edges[edge]["prob"] = self.prob_class(
                reduce(
                    multiply, [data["prob"] for data in edges_for_same_nodes.values()]
                ),
                level=self.poly_level,
            )
            graph.edges[edge]["multiplicity"] = len(multy_graph[edge[0]][edge[1]])
        return graph

    def total_weight(self):
        """Return the total weight of the graph nodes"""
        return np.sum(
            [self.graph.nodes[node]["weight"] for node in self.graph.nodes], axis=0
        )

    def edge_total_weight(self):
        """Return the total reliability weight stored on graph edges."""
        return np.sum(
            [data.get("edge_weight", 0.0) for _, _, data in self.graph.edges(data=True)],
            axis=0,
        )

    def decompose(self, *, include_generalized_chains: bool = True):
        """Return an edge-connectivity decomposition of this reliability graph.

        The decomposition is intentionally limited to simple input graphs. It
        rejects user-supplied MultiGraphs and also rejects cases where source
        contraction creates parallel edges.
        """
        from .rel_decomposition import decompose_graph_rel

        return decompose_graph_rel(
            self,
            include_generalized_chains=include_generalized_chains,
        )

    def _raise_if_edge_weights_for_exact_methods(self):
        if float(self.edge_total_weight()) != 0.0:
            raise NotImplementedError(
                "edge_weight is currently supported only by calc_rel_simulation()"
            )


def max_key(graph: nx.MultiGraph, u, v):
    """Find the maximum key value for edges between two nodes.

    ###Args:
        graph (nx.MultiGraph): The graph being analyzed.
        u: The first node.
        v: The second node.

    ###Returns:
        int: The maximum key value for edges between the two nodes.
    """
    if graph.has_edge(u, v):
        return max(graph[u][v])
    return 0


def contracted_nodes_list(graph: nx.MultiGraph, nodes: ArrayLike) -> nx.MultiGraph:
    """Contract a list of nodes into a single node in the graph, preserving self-loops.

    ###Args:
        graph (nx.MultiGraph): The graph being analyzed.
        nodes (ArrayLike): The list of nodes to contract.

    ###Returns:
        nx.MultiGraph: The graph with the nodes contracted.
    """
    res_graph = nx.MultiGraph(graph)
    for node in nodes[1:]:
        nx.contracted_nodes(res_graph, nodes[0], node, self_loops=True, copy=False)
    return res_graph


def non_self_edge(graph: nx.MultiGraph, v) -> tuple:
    """Find a non-self edge connected to the specified node.

    ###Args:
        graph (nx.MultiGraph): The graph being analyzed.
        v: The node to find a non-self edge for.

    ###Returns:
        tuple: A non-self edge connected to the specified node.
    """
    for edge in graph.edges(v):
        if edge[0] != edge[1]:
            return edge
    return None


def rename_nodes_prob_to_original_labels(
    gr: GraphRel, nodes_probs: Dict[int, PROBS.Prob]
) -> Dict:
    """Rename the nodes prob from gr.graph edges to gr.original_graph edges"""
    for node in gr.original_graph.nodes:
        node_int_label = gr.original_graph.nodes[node]["node_int_label"]
        if node_int_label in nodes_probs:
            gr.original_graph.nodes[node]["rel"] = nodes_probs[node_int_label]
        else:
            gr.original_graph.nodes[node]["rel"] = gr.prob_class([0.0])
    return nx.get_node_attributes(gr.original_graph, "rel")


##### init #####


def define_prob_class_and_init_edges_prob(
    edges_prob: Dict[tuple, PROBS.Prob], max_fail: int
) -> Tuple[ABCMeta, Dict[tuple, PROBS.Prob]]:
    """Get the edges_prob of the graph edges and thier probability class.
    Return the edge prob according to the prob_class. The main use of the function
    is in Poly prob_class. The degree of the polynom is set to ensure a fixes size array
    in the reliability calculations and speed up the computations.

    Parameters
    ----------
    edges_prob : Dict[tuple, PROBS.Prob]
        The failure probability of each edge
    max_fail : int
        The maximal number of failing edges(used to approximate the reliabilty)

    Returns
    -------
    Tuple[ABCMeta, Dict[tuple, PROBS.Prob]]
        A tuple of the prob_class ofthe edges and the new edges_prob
    """
    _validate_edges_prob(edges_prob)
    # If the prob is a Poly, define the degree and the maximal failing of the polynom
    prob_class = next(iter(edges_prob.values())).__class__
    if prob_class == PROBS.Poly:
        degree = get_degree(list(edges_prob.values()))
        if max_fail is not None:
            degree = min([degree, max_fail])
        prob_class = PROBS.create_poly_class(degree=degree)
        new_edges_prob = {edge: prob_class(prob) for edge, prob in edges_prob.items()}
        return prob_class, new_edges_prob
    return prob_class, edges_prob


def _validate_edges_prob(edges_prob: Dict[tuple, PROBS.Prob]):
    """Check that all the edges are from the same prob class"""
    # validate that all the probs are instance of Prob
    if not all([isinstance(prob, PROBS.Prob) for prob in edges_prob.values()]):
        raise TypeError("All the probs on edges_prob must be from type Prob")
    # validate that all the probs are from the same class
    prob_class_names = set(prob.__class__.__name__ for prob in edges_prob.values())
    if len(prob_class_names) > 1:
        raise TypeError("All the edges_prob must be from the same class")


def saidi_from_nodes_prob_and_weight(nodes_prob: dict, nodes_weight: dict) -> float:
    """
    Get the discinnection probability of each node and the weight of each node.
    Return the SAIDI index of the graph be summing saidi = sum(nodes_prob[v] * nodes_weight[v]) / total_weight
    """
    disconnection_weight = np.sum(
        [nodes_prob[node] * nodes_weight[node] for node in nodes_prob], axis=0
    )
    total_weight = np.sum(list(nodes_weight.values()), axis=0)
    return disconnection_weight / total_weight


def edges_rel_in_original_graph(gr: GraphRel, edges_prob: dict) -> dict:
    """Transform the edges reliability from gr.graph to gr.original_graph edges"""
    return {
        gr.graph.edges[edge]["original_edge"]: prob for edge, prob in edges_prob.items()
    }
