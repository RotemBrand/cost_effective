"""
The tree_decomposition.py file is responsible to operate the tree decomposition-based
algorithm to calculate the graph reliability indexes.
"""
from abc import ABC
import numpy as np
from typing import List
from numpy.typing import ArrayLike
import networkx as nx
import indexes.partition_collection as PC
from indexes.probs import Poly
from functools import reduce
import indexes.utilities as utilities
import networkx.algorithms.approximation.treewidth as treewidth


### atr calc #####

def calc_atr(gr: "GR.GraphRel") -> Poly:
    bugs_order = gr.td.bugs_order
    for bug in bugs_order:
        in_edges = [edge for edge in gr.td.in_edges(bug) if edge[0] != edge[1]]
        out_edge = next(iter((gr.td.out_edges(bug))))
        if out_edge[0] != out_edge[1]:
            base_partition_collection = (
                gr.td.edges[out_edge]["in_partitions_collection"]
                * gr.td.nodes[bug]["boundary_partitions_collection"]
            )
        else:
            base_partition_collection = gr.td.nodes[bug][
                "boundary_partitions_collection"
            ]
        # multiply
        for in_edge in in_edges:
            partitions_for_product = [base_partition_collection] + [
                gr.td.edges[in_edge2]["out_partitions_collection"]
                for in_edge2 in in_edges
                if in_edge2 != in_edge
            ]
            in_partitions_collection = PC.cartesian_product_partitions_collection_list(
                partitions_for_product
            )
            in_partitions_collection = in_partitions_collection.project(in_edge[0])
            gr.td.edges[in_edge]["in_partitions_collection"] = in_partitions_collection


### saidi main algo steps ####


def set_the_final_partitions_collection_for_each_bug(gr, run_out_pc: bool):
    r"""run the tree decomposition algorithm to calculate the saidi\nodes_rel\edges_rel
    indexes. The output of the algorithm is to set for each but of the tree decomposition
    the attribute "final_partitions_collection" that represent for each bug the probability for each partition
    A partition of the bug nodes represent the conncected componnnets of the bug nodes.

    Parameters
    ----------
    gr : GR.GraphRel
        A GraphRel object
    """
    calculate_bugs_partitions_collection(gr)
    if run_out_pc:
        calculate_out_partitions_collection(gr)
    calculate_in_partitions_collections(gr)
    calculate_bug_final_partitions_collection(gr)


def calculate_bugs_partitions_collection(gr: "GR.GraphRel"):
    """Calculate for each bug the partitions_collection attribute which is
    the probability for each connected components of the bug subgraph.
    The calculations uses the pivot method

    Parameters
    ----------
    gr : GR.GraphRel
        A Graph Rel object
    """
    for bug, bug_data in gr.td.nodes(data=True):
        gr.td.nodes[bug]["partitions_collection"] = get_partitions_probability_of_graph(
            graph=bug_data["subgraph"],
            # source_nodes=bug,
            source_nodes=bug_data["boundary"],
            prob_class=gr.prob_class
        )
        gr.td.nodes[bug]["boundary_partitions_collection"] = gr.td.nodes[bug][
            "partitions_collection"
        ].project(bug_data["boundary"])


def calculate_out_partitions_collection(gr: "GR.GraphRel"):
    """Get for each edge of the tree decomposition the out_partitions_collection attribute
    the out_partitions_collection is the probability for each connected componnents of the
    edge intersection nodes on the graph that before the edge

    Parameters
    ----------
    gr : _type_
        _description_
    """
    bugs_order = gr.td.bugs_order[::-1]
    for bug in bugs_order:
        in_edges = [edge for edge in gr.td.in_edges(bug) if edge[0] != edge[1]]
        out_edge = next(iter((gr.td.out_edges(bug))))
        # multiply
        in_partitions_collection_list = [
            gr.td.edges[edge]["out_partitions_collection"] for edge in in_edges
        ]
        in_partitions_collection_list.append(
            gr.td.nodes[bug]["boundary_partitions_collection"]
        )
        out_partitions_collection = PC.cartesian_product_partitions_collection_list(
            in_partitions_collection_list
        )
        out_partitions_collection = out_partitions_collection.project(out_edge[1])
        gr.td.edges[out_edge]["out_partitions_collection"] = out_partitions_collection


def calculate_in_partitions_collections(gr: "GR.GraphRel"):
    """
    Set the in_partitions_collections attribute for each edge of the td graph
    """
    bugs_order = gr.td.bugs_order
    for bug in bugs_order:
        in_edges = [edge for edge in gr.td.in_edges(bug) if edge[0] != edge[1]]
        out_edge = next(iter((gr.td.out_edges(bug))))
        if out_edge[0] != out_edge[1]:
            base_partition_collection = (
                gr.td.edges[out_edge]["in_partitions_collection"]
                * gr.td.nodes[bug]["boundary_partitions_collection"]
            )
        else:
            base_partition_collection = gr.td.nodes[bug][
                "boundary_partitions_collection"
            ]
        # multiply
        for in_edge in in_edges:
            partitions_for_product = [base_partition_collection] + [
                gr.td.edges[in_edge2]["out_partitions_collection"]
                for in_edge2 in in_edges
                if in_edge2 != in_edge
            ]
            in_partitions_collection = PC.cartesian_product_partitions_collection_list(
                partitions_for_product
            )
            in_partitions_collection = in_partitions_collection.project(in_edge[0])
            gr.td.edges[in_edge]["in_partitions_collection"] = in_partitions_collection


def calculate_bug_final_partitions_collection(gr: "GR.GraphRel"):
    """
    Set the bug final_partitions_collection attribute for each node of gr.td
    """
    for bug in gr.td.nodes:
        in_edges = [edge for edge in gr.td.in_edges(bug) if edge[0] != edge[1]]
        out_edge = next(iter((gr.td.out_edges(bug))))
        out_pc = (
            [gr.td.edges[out_edge]["in_partitions_collection"]]
            if out_edge[0] != out_edge[1]
            else []
        )
        partitions_collections = (
            [gr.td.nodes[bug]["partitions_collection"]]
            + out_pc
            + [
                gr.td.edges[in_edge]["out_partitions_collection"]
                for in_edge in in_edges
                if in_edge[0] != in_edge[1]
            ]
        )
        gr.td.nodes[bug]["final_partitions_collection"] = (
            PC.cartesian_product_partitions_collection_list(partitions_collections)
        )


### get the nodes and edges disconnection prob from the final_partitions_collection ###
def get_nodes_disconnection_prob_from_pc(
    gr: "GR.GraphRel", pc: PC.PartitionsCollection, source
) -> dict:
    """
    Get a gr object after apply the set_the_final_partitions_collection_for_each_bug function
    with the final_partitions_collection attribute for each bug. Get the disconnection probability of each
    of the original graph nodes
    """
    pc_non_fset = pc.change_is_fset(False)
    connection_probs = {}
    for partition, prob in pc_non_fset.items():
        source_comp = frozenset()
        # find the source component of the partition
        for comp in partition:
            if source in comp:
                source_comp = comp
                break
        # update the connection prob of each node
        for node in source_comp:
            if node not in connection_probs:
                connection_probs[node] = prob
            else:
                connection_probs[node] += prob
    # connection prob to disconnection prob
    connection_probs[source] = gr.prob_class(1)
    disconnection_prob = {
        node: gr.prob_class(1) - prob for node, prob in connection_probs.items()
    }
    return disconnection_prob


def key_edge(edge: tuple) -> tuple:
    """Get a unique key foe each edge by sort its nodes"""
    return tuple(sorted(edge[:2]) + [edge[2]])


def get_edges_connections_prob_from_pc(
    pc: PC.PartitionsCollection, source, edges: list
) -> dict:
    """
    Get a gr object after apply the set_the_final_partitions_collection_for_each_bug function
    with the final_partitions_collection attribute for each bug. Get the reaching probability of each
    of the original graph edges
    """
    # init
    pc_non_fset = pc.change_is_fset(False)
    connection_probs = {
        key_edge(edge): {
            frozenset(subset): pc.prob_class([0]) for subset in utilities.all_subsets(edge[:2])
        }
        for edge in edges
    }
    for partition, prob in pc_non_fset.items():
        # find the source component
        source_comp = frozenset([source])
        for comp in partition:
            if source in comp:
                source_comp = comp
                break
        source_comp = frozenset(source_comp)
        # add the reaching probabilities of each edge
        for edge in edges:
            intersection = frozenset(edge[:2]).intersection(source_comp)
            connection_probs[key_edge(edge)][intersection] += prob
    return connection_probs


def get_nodes_disconnection_prob_from_td(gr: "GR.GraphRel") -> dict:
    """
    Get a gr object after apply the set_the_final_partitions_collection_for_each_bug function
    with the final_partitions_collection attribute for each bug. Get the disconnection probability of each
    of the original graph nodes
    """
    disconnection_prob = {}
    for bug, bug_data in gr.td.nodes.items():
        bug_disconnection_probs = get_nodes_disconnection_prob_from_pc(
            gr, bug_data["final_partitions_collection"], source=gr.source
        )
        for node, prob in bug_disconnection_probs.items():
            if node not in disconnection_prob:
                gr.graph.nodes[node]["disconnection_prob"] = prob
                disconnection_prob[node] = prob
    return disconnection_prob


def set_edges_disconnection_prob(gr: "GR.GraphRel") -> dict:
    """"""
    edges_connections_prob = {}
    for bug, bug_data in gr.td.nodes.items():
        bug_disconnection_probs = get_edges_connections_prob_from_pc(
            bug_data["final_partitions_collection"],
            source=gr.source,
            edges=list(bug_data["subgraph"].edges),
        )
        for edge, prob_dict in bug_disconnection_probs.items():
            if edge not in edges_connections_prob:
                edges_connections_prob[edge] = get_edge_rel_without_the_edge(
                    prob_dict=prob_dict,
                    edge=edge,
                    edge_fail_prob=gr.graph.edges[edge]["prob"],
                )
    nx.set_edge_attributes(gr.graph, edges_connections_prob, "connection_probs")

    return edges_connections_prob


def get_edge_rel_without_the_edge(prob_dict: dict, edge: tuple, edge_fail_prob) -> dict:
    """
    Get the reaching probabilities of edges in the graph G and
    return the edge reaching probabilities in the graph G/e for each edge e
    """
    # hande self edges
    if edge[0] == edge[1]:
        return prob_dict.copy()
    new_prob_dict = {}
    new_prob_dict[frozenset()] = prob_dict[frozenset()]
    new_prob_dict[frozenset(edge[0:1])] = (
        prob_dict[frozenset(edge[0:1])] // edge_fail_prob
    )
    new_prob_dict[frozenset(edge[1:2])] = (
        prob_dict[frozenset(edge[1:2])] // edge_fail_prob
    )
    new_prob_dict[frozenset(edge[:2])] = (
        -np.sum(list(new_prob_dict.values()), axis=0) + 1
    )
    return new_prob_dict


#### partitions collection from graph ####


def get_partitions_probability_of_graph(
    graph: nx.Graph, source_nodes: ArrayLike, prob_class: ABC
) -> PC.PartitionsCollection:
    """
    Get a graph and calculate the probability of each connected components of the graph
    using the pivot method
    """
    nx.set_node_attributes(
        graph, {node: frozenset([node]) for node in graph.nodes}, "component"
    )
    partitions_collection = get_partitions_prob_recursive(
        graph.copy(), source_nodes, prob_class
    )
    return partitions_collection.change_is_fset(True)


def get_partitions_prob_recursive(
    graph: nx.MultiGraph, projection_nodes: frozenset, prob_class: ABC
) -> PC.PartitionsCollection:
    """
    The recursive function of the pivot method used to calculate the probability for each
    connected components of the graph
    """
    edge = get_non_self_edge_of_projection_nodes(graph, projection_nodes)
    # edge = get_non_self_edge(graph)
    if edge is None:
        components = frozenset(
            node_data["component"]
            for node_data in graph.nodes.values()
            if len(node_data["component"]) > 1
        )
        return PC.PartitionsCollection(
            {components: prob_class(1)}, prob_class=prob_class, is_fset=False
        )
    edge = edge[:2]
    edge_data = graph.get_edge_data(edge[0], edge[1])
    graph_contracted_edge = contracted_edge(graph, edge)
    graph_delete_edge = graph.copy()
    graph_delete_edge.remove_edges_from(
        [list(edge) + [key] for key in edge_data.keys()]
    )

    partitions_collection_contracted = get_partitions_prob_recursive(
        graph_contracted_edge, projection_nodes, prob_class
    )
    partitions_collection_deleted = get_partitions_prob_recursive(
        graph_delete_edge, projection_nodes, prob_class
    )
    failing_prob = np.prod([data["prob"] for data in edge_data.values()], axis=0)
    return partitions_collection_contracted.multiply_prob(
        prob_class(1) - failing_prob
    ) + partitions_collection_deleted.multiply_prob(failing_prob)


def get_non_self_edge(graph: nx.MultiGraph) -> tuple:
    """Output a non-self edge of the graph"""
    for edge in graph.edges:
        if edge[0] != edge[1]:
            return edge
    return None

def get_non_self_edge_of_projection_nodes(graph: nx.MultiGraph, projection_nodes: frozenset) -> tuple:
    res_edge = None
    for node in projection_nodes:
        if node in graph.nodes:
            for edge in graph.edges(node):
                if edge[0] != edge[1]:
                    res_edge = edge
                    break
    return res_edge

def contracted_edge(graph, edge):
    """
    contract the graph with the edge nodes
    """
    (u, v) = edge
    res = nx.contracted_edge(graph, edge, self_loops=False)
    res.nodes[u]["component"] = graph.nodes[u]["component"].union(
        graph.nodes[v]["component"]
    )
    return res


##### create the tree decomposition graph #####


def get_td_with_attributes(gr: "GR.GraphRel") -> nx.DiGraph:
    """Create the tree decomposition used for the td calculations
    and prepare it for computations be set different edge and node attributes:

    node_attributes:
    -------------
        - boundary(frozenset): the graphs that are on the bug boundery
        - subgraph(frozenset): the subgraph of the original graph that the bug represent
        - partitions_collection(PartitionCollection): The partition collection that represent
        all the ways that the subgraph nodes are connected and their probability

    edges_attributes:
    -------------
        - intersection(frozenset): the nodes intersection of the two bugs connected by the edge
    Parameters
    ----------
    gr : GraphRel
        A GraphRel object

    Returns
    -------
    _type_
        A graph that represent the tree decomposition of the orignal graph
    """
    td = treewidth.treewidth_min_fill_in(nx.Graph(gr.graph))[1]
    # td = get_td(gr.graph, gr.source)
    # add_bugs_order_delete(td, gr.source)
    add_source_node_to_each_bug(td, gr.source)
    td = tree_to_digraph(td, next(iter(td.nodes)))
    add_intersecion_nodes_for_td_edges(td)
    add_boundary_nodes_to_bugs(td)
    add_subgraph_to_bugs(td, gr.graph)
    return td


def get_td(graph: nx.Graph, source):
    print("Warning: use skeleton td")
    skeleton = utilities.get_skeleton_graph(graph, [source])
    td_skeleton: nx.Graph = treewidth.treewidth_min_fill_in(nx.Graph(skeleton))[1]
    # get the bug of each chain
    chain_to_bug = {}
    for bug in td_skeleton.nodes:
        subgraph = nx.subgraph(skeleton, bug)
        for chain in subgraph.edges:
            if chain not in chain_to_bug:
                chain_to_bug[chain] = bug
    # add the chains subgraphs to the td
    for chain, bug in chain_to_bug.items():
        chain_subgraph = skeleton.edges[chain]['subgraph']
        nodes_order = list(nx.dfs_preorder_nodes(chain_subgraph, chain[0]))
        sublists = [frozenset(nodes_order[i:i + 3]) for i in range(0, len(nodes_order), 2)]
        sublists = [sublist.union(frozenset([nodes_order[0], nodes_order[-1]])) for sublist in sublists]
        td_skeleton.add_nodes_from(sublists)
        td_skeleton.add_edges_from([(sublists[i], sublists[i + 1]) for i in range(len(sublists) - 1)])
        td_skeleton.add_edge(sublists[0], bug)
    return td_skeleton

        

def add_source_node_to_each_bug(td: nx.Graph, source) -> None:
    """add the source node to each bug of the td"""
    new_bugs = {bug: frozenset.union(bug, {source}) for bug in td.nodes}
    nx.relabel_nodes(td, new_bugs, copy=False)


def tree_to_digraph(tree: nx.Graph, source) -> nx.DiGraph:
    """Transform a tree to a digraph by tunning bfs starting from the tree source"""
    # Create an empty directed graph
    digraph = nx.DiGraph()
    # Start from the source node and do a breadth-first traversal
    for parent, child in nx.bfs_edges(tree, source):
        digraph.add_edge(parent, child)
    bugs_order = list(nx.topological_sort(digraph))
    digraph.add_edge(source, source)
    digraph = digraph.reverse()
    digraph.bugs_order = bugs_order
    return digraph


def add_subgraph_to_bugs(td: nx.DiGraph, graph: nx.Graph):
    """add the subgraph attributes to each of the td bugs"""
    for bug in td.nodes:
        td.nodes[bug]["subgraph"] = nx.MultiGraph(nx.subgraph(graph, bug))
    delete_edges_from_multiple_td_nodes(td)


def delete_edges_from_multiple_td_nodes(td: nx.DiGraph):
    """
    For each edge that exist in multiple bugs subgraph,
    remove the edge from all the subgraph except from one of them
    """
    existing_edges = set()
    for bug in td.bugs_order[::-1]:
        subgraph = td.nodes[bug]["subgraph"]
        edge_to_remove = []
        for edge in subgraph.edges:
            if frozenset(edge) in existing_edges:
                edge_to_remove.append(edge)
            else:
                existing_edges.add(frozenset(edge))
        subgraph.remove_edges_from(edge_to_remove)


def add_intersecion_nodes_for_td_edges(td: nx.DiGraph):
    """add the intersection attribute to the td edges"""
    for edge in td.edges:
        td.edges[edge]["intersection"] = edge[0].intersection(edge[1])


def add_boundary_nodes_to_bugs(td: nx.DiGraph):
    """Add for each bug in the td, the intersection of the bug with its neighboors"""
    for bug in td.nodes:
        neighboors = list(td.predecessors(bug)) + list(td.successors(bug))
        neighboors_union = reduce(frozenset.union, neighboors)
        td.nodes[bug]["boundary"] = neighboors_union.intersection(bug)


def add_bugs_order(td: nx.DiGraph):
    """Add order to the td bugs by running topological_sort"""
    td.bugs_order = list(nx.topological_sort(td))[::-1]
