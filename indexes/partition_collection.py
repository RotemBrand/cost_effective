"""
Make partitions collections computations

Partition:
--------
    * A partition of a set represent a possible connected components of the set

FSet:
-------
    * Fset is a set that it represented by integers. The integer represent all the elements i that has one in the i bit

PartitionsCollection:
------------------
    * Represent the distribution of all the possible partitions of a set
"""

import numpy as np
from collections import UserDict, deque

# import indexes.graph_rel as gr
from numpy.typing import ArrayLike
from typing import Dict, Union, FrozenSet, List
import networkx as nx
from itertools import combinations
from copy import copy
import indexes.probs as PROBS


class FSet(int):
    """
    set is a set that it represented by integers. The integer represent all the elements i that has one in the i bit
    """

    @classmethod
    def from_set(cls, s: "FSet") -> "FSet":
        """## Get a set and return a number that represent the set.
        The number i is in the set if the i'th bit in the number is 1

        ### Args:
            - `s (FrozenSet)`: A set of numbers

        ### Returns:
            - `int`: The number that represent the set
        """
        fset = 0
        for x in s:
            fset |= 1 << x
        return fset

    @classmethod
    def to_set(cls, fset: "FSet") -> frozenset:
        """
        Gets a nunber that represent a set: The number i is in the set if the i'th bit in the number is 1
        Return the represented set
        Parameters
        ----------
        fset : FSet

        Returns
        -------
        frozenset
            A set that conatain the index of all the 1 bits of the number

        """
        mask = 1
        i = 0
        res_list = []
        while mask <= fset:
            if fset & mask:  # check if the i'th bit is in num
                res_list.append(i)
            mask <<= 1
            i += 1
        return frozenset(res_list)

    @classmethod
    def is_singleton_or_empty(cls, fset: "FSet") -> bool:
        """
        Parameters
        ----------
        fset : FSet

        Returns
        -------
        bool
            Check if the fset contain only one or zero bit of 1. It also check
            if the number represent a singelton

        """
        return fset & (fset - 1) == 0

    @classmethod
    def sets_to_fsets(cls, set_of_sets: FrozenSet[FrozenSet[int]]) -> FrozenSet["FSet"]:
        """
        cast a set of sets on number to set of fast sets. Each set in the big set represnt by a number.
        Parameters
        ----------
        set_of_sets : set of sets

        Returns
        -------
        frozenset of numbers
            For each set in the input return the number represent the set

        """
        return frozenset(FSet.from_set(s) for s in set_of_sets)

    @classmethod
    def fsets_to_sets(cls, fsets: FrozenSet["FSet"]) -> FrozenSet[FrozenSet[int]]:
        """
        cast a set of numbers on set of sets. Each number in the set represent a set.
        Parameters
        ----------
        fsets : set of FSets(numbers)

        Returns
        -------
        frozen set of frozensets
            For each number in the input set return the matching set.

        """
        return frozenset(FSet.to_set(num) for num in fsets)

    @classmethod
    def intersects(cls, num1: "FSet", num2: "FSet") -> bool:
        """
        Return True if the sets that represented by numbers are intersected
        We do it by making bitwise and
        """
        return num1 & num2 > 0

    @classmethod
    def intersection(cls, num1: "FSet", num2: "FSet") -> bool:
        """
        Return the intersection of two fsets
        """
        return num1 & num2

    @classmethod
    def union(cls, num1: "FSet", num2: "FSet") -> "FSet":
        """
        Return the union of the fast sets
        """
        return num1 | num2

    @classmethod
    def to_minimal_form(cls, partition: FrozenSet["FSet"]) -> FrozenSet["FSet"]:
        components_queue = deque()
        for fset in partition:
            counter = 0
            max_counter = len(components_queue)
            while counter < max_counter:
                fset_res = components_queue.popleft()
                if FSet.intersection(fset_res, fset):
                    fset = FSet.union(fset_res, fset)
                else:
                    components_queue.append(fset_res)
                counter += 1
            components_queue.append(fset)
        return frozenset(components_queue)

    @classmethod
    def remove_singletons(cls, fset_partition: FrozenSet["FSet"]) -> FrozenSet["FSet"]:
        """
        Remove singeltons from a partition of fsets

        parameters
        ----------
        partition : set of Fsets

        Returns
        -------
        frozenset
            A partition without singeltons(sets from size one). It important
            because nodes in a singelton are disconnected from any node
        """
        return frozenset(
            [
                component
                for component in fset_partition
                if not FSet.is_singleton_or_empty(component)
            ]
        )


def remove_singleton(partition: FrozenSet[FrozenSet[int]]) -> FrozenSet[FrozenSet[int]]:
    """
    Remove singeltons from a partition

    parameters
    ----------
    partition : set of sets

    Returns
    -------
    frozenset
        A partition without singeltons(sets from size one). It important
        because nodes in a singelton are disconnected from any node
    """
    return frozenset(
        [frozenset(component) for component in partition if len(component) > 1]
    )


def add_dict_to_dict(
    d_target: dict, d_update: dict, inplace: bool = True
) -> Union[None, dict]:
    """
    Adds the values of one dictionary (`d_update`) to another (`d_target`),
    either modifying the target dictionary in place or returning a new dictionary.

    If a key in `d_update` exists in `d_target`, their values are combined using addition (`+`).
    If a key does not exist in `d_target`, it is copied from `d_update`.

    Parameters
    ----------
    d_target : dict
        The target dictionary to which values from `d_update` will be added.
    d_update : dict
        The dictionary containing values to add to `d_target`.
        Must have an `items` method (like a standard dictionary).
    inplace : bool, optional
        If True, modifies `d_target` directly.
        If False, returns a new dictionary with the combined values, leaving `d_target` unchanged.
        Default is True.

    Returns
    -------
    Union[None, dict]
        If `inplace` is True, returns None (modifies `d_target` in place).
        If `inplace` is False, returns a new dictionary with the combined values.

    Raises
    ------
    TypeError
        If `d_update` does not have an `items` method.
    """

    if not hasattr(d_update, "items"):
        raise TypeError("d_updates must have .items function")
    d_res = d_target if inplace else copy(d_target)
    for key, value in d_update.items():
        if key in d_res:
            d_res[key] += value
        else:
            d_res[key] = copy(value)
    if not inplace:
        return d_res


COUNTER = 0  # TODO: delete


def cartesian_product_partition(
    partition1: FrozenSet[FSet], partition2: FrozenSet[FSet]
) -> FrozenSet[FSet]:
    """
    Computes the Cartesian product of two partitions represented as `FSet` objects.

    Parameters
    ----------
    partition1 : FrozenSet[FSet]
        The first partition, where each component is an `FSet`.
    partition2 : FrozenSet[FSet]
        The second partition, where each component is an `FSet`.

    Returns
    -------
    FrozenSet[FSet]
        A new partition representing the Cartesian product of `partition1` and `partition2`.
        Overlapping components between the partitions are merged.

    Notes
    -----
    - The function iteratively merges overlapping components between the two partitions.
    - The input partitions are processed in ascending order of size for efficiency.
    - Uses `FSet.union` to merge overlapping components and `FSet.intersects` to detect overlaps.
    - The global `COUNTER` variable is incremented for every component comparison for debugging purposes.

    Examples
    --------
    Given two partitions:
    - `partition1 = frozenset({FSet.from_set({0, 1}), })`
    - `partition2 = frozenset({FSet.from_set({1, 2}), FSet.from_set({3, 4})})`

    The result will be:
    - `frozenset({FSet.from_set({0, 1, 2}), FSet.from_set({3, 4})})`
    """
    global COUNTER
    partition_small, partition_large = sorted([partition1, partition2], key=len)
    components_queue = deque(partition_large)
    for fset2 in partition_small:
        counter = 0
        max_counter = len(components_queue)
        while counter < max_counter:
            COUNTER += 1
            fset_res = components_queue.popleft()
            if FSet.intersects(fset_res, fset2):
                fset2 = FSet.union(fset_res, fset2)
            else:
                components_queue.append(fset_res)
            counter += 1
        components_queue.append(fset2)
    return frozenset(components_queue)


class PartitionsCollection:
    """
    Represents a collection of partitions with associated probabilities.

    Each partition is either represented as sets or as `FSet` objects (bitwise representations),
    depending on the `is_fset` attribute. Provides operations for transformation, combination,
    and manipulation of these partitions and their probabilities.

    Attributes
    ----------
    _pc_dict : dict
        The dictionary where keys are partitions and values are their associated probabilities.
    prob_class : PROBS.Prob
        The probability class used for computations (e.g., addition, multiplication).
    _is_fset : bool
        Indicates whether partitions are represented as `FSet` objects (True) or regular sets (False).
    """

    def __init__(
        self,
        pc_dict: dict,
        prob_class: PROBS.Prob,
        is_fset: bool,
        is_zeta: bool=False,
        remove_singletons: bool = False,
    ):
        """
        Initialize a PartitionsCollection instance.

        Parameters
        ----------
        pc_dict : dict
            Dictionary where keys are partitions and values are associated probabilities.
        prob_class : PROBS.Prob
            Probability class to manage probability operations.
        is_fset : bool
            If True, partitions are represented as `FSet` objects; otherwise, as sets.
        is_zeta: bool
            If True, the partition represent a zeta transform
        remove_singletons : bool, optional
            If True, removes singleton partitions (default is False).
        """
        # TODO: validate dict, validate singeltons
        self._pc_dict = pc_dict
        self.prob_class = prob_class
        self._is_fset = is_fset
        self._is_zeta = is_zeta
        if remove_singletons:
            self.remove_singletons()

    @property
    def is_fset(self) -> bool:
        """
        Indicates whether partitions are represented as `FSet` objects.

        Returns
        -------
        bool
            True if partitions are `FSet` objects, False otherwise.
        """
        return self._is_fset

    @is_fset.setter
    def is_fset(self, is_fset_option: bool):
        """
        Sets whether partitions should be represented as `FSet` objects.

        Parameters
        ----------
        is_fset_option : bool
            True to represent partitions as `FSet`, False for sets.

        Notes
        -----
        Updates the internal partition dictionary representation accordingly.
        """
        if is_fset_option != self._is_fset:
            if self._is_fset:
                self._pc_dict = {
                    FSet.fsets_to_sets(partition): self.prob_class(prob)
                    for partition, prob in self._pc_dict.items()
                }
            else:
                self._pc_dict = {
                    FSet.sets_to_fsets(partition): copy(prob.prob)
                    for partition, prob in self._pc_dict.items()
                }
            self._is_fset = not self._is_fset

    @property
    def pc_dict(self) -> dict:
        """
        Returns the internal dictionary of partitions and probabilities.

        Returns
        -------
        dict
            Dictionary of partitions and their associated probabilities.
        """
        return self._pc_dict

    def change_is_fset(self, is_fset: bool, inplace: bool = False):
        """
        Converts the representation of partitions between sets and `FSet`.

        Parameters
        ----------
        is_fset : bool
            True to convert to `FSet` representation, False for sets.
        inplace : bool, optional
            If True, modifies the current object; otherwise, returns a new object.

        Returns
        -------
        PartitionsCollection or None
            A new PartitionsCollection if `inplace` is False; otherwise, None.
        """

        if inplace:
            self.is_fset = is_fset
        else:
            new_pc = PartitionsCollection(
                pc_dict=self._pc_dict.copy(),
                prob_class=self.prob_class,
                is_fset=self.is_fset,
            )
            new_pc.is_fset = is_fset
            return new_pc

    def remove_singletons(self):
        """
        Removes singleton partitions from the collection.

        Notes
        -----
        The function modifies the internal partition dictionary by removing partitions
        that contain single components.
        """
        remove_singletons_func = (
            FSet.remove_singletons if self._is_fset else remove_singleton
        )
        new_pc_dict = {}
        for partition, prob in self._pc_dict.items():
            new_partition = remove_singletons_func(partition)
            add_dict_to_dict(new_pc_dict, {new_partition: prob})
        self._pc_dict = new_pc_dict

    def copy(self) -> "PartitionsCollection":
        """
        Creates a deep copy of the current PartitionsCollection.

        Returns
        -------
        PartitionsCollection
            A new PartitionsCollection object with the same data.
        """
        return PartitionsCollection(
            pc_dict=self.pc_dict.copy(),
            prob_class=self.prob_class,
            is_fset=self.is_fset,
        )

    def __len__(self):
        """
        Returns the number of partitions in the collection.

        Returns
        -------
        int
            Number of partitions.
        """
        return len(self._pc_dict)

    def __str__(self):
        items = []
        for partition, prob in self._pc_dict.items():
            if self.is_fset:
                partition_copy = tuple(sorted([comp] for comp in partition))
            else:
                partition_copy = tuple(sorted(sorted(comp) for comp in partition))
            items.append([partition_copy, prob])
        str_list = [f"\t{partition}: {prob}" for partition, prob in items]
        return "{\n" + "\n".join(str_list) + "\n}\n" + f"len = {len(self)}"

    def __repr__(self):
        return self.__str__()

    def get_all_elements(self) -> List[int]:
        """
        Retrieves all unique elements present in the partitions.

        Returns
        -------
        List[int]
            List of unique elements in the partitions.

        Raises
        ------
        ValueError
            If the collection is using `FSet` representation.
        """

        if self.is_fset:
            raise ValueError(
                "get_all_elements is not valid for fset partition collection"
            )
        return set(
            [elm for partition in self.pc_dict for comp in partition for elm in comp]
        )

    def project(self, projection_set: ArrayLike) -> "PartitionsCollection":
        """
        Projects the partitions onto a specified subset of elements.

        Parameters
        ----------
        projection_set : ArrayLike
            The set of elements to project the partitions onto.

        Returns
        -------
        PartitionsCollection
            A new PartitionsCollection containing partitions restricted to the `projection_set`.

        Notes
        -----
        - Converts the partitions to `FSet` representation if not already in that form.
        - Singleton partitions are removed from the result.
        """
        # set the pc to be fset
        old_is_fset = self._is_fset
        self.is_fset = True
        projection_fset = FSet.from_set(projection_set)
        # make projection
        new_pc_dict = {}
        for partition_fset, prob in self._pc_dict.items():
            new_partition_fset = FSet.remove_singletons(
                FSet.intersection(fcomp, projection_fset) for fcomp in partition_fset
            )
            add_dict_to_dict(new_pc_dict, {new_partition_fset: prob}, inplace=True)
        # return the new pc
        new_pc = PartitionsCollection(
            pc_dict=new_pc_dict, prob_class=self.prob_class, is_fset=True
        ).change_is_fset(old_is_fset)
        return new_pc

    def __mul__(self, other: "PartitionsCollection") -> "PartitionsCollection":
        """
        Multiplies two PartitionsCollection objects, creating a Cartesian product of their partitions.

        Parameters
        ----------
        other : PartitionsCollection
            The other PartitionsCollection to multiply with the current collection.

        Returns
        -------
        PartitionsCollection
            A new PartitionsCollection representing the Cartesian product of the two inputs.

        Raises
        ------
        ValueError
            If either of the collections has `is_fset` set to `False`.

        Notes
        -----
        - Both collections must be in `FSet` format (`is_fset=True`) to perform multiplication.
        - Updates probabilities using `Prob.fast_prod`.
        - Filters out partitions with zero probabilities if `Prob.need_to_check_zero_events` is `True`.
        - Temporarily changes the `is_fset` state of the input collections if necessary, reverting them afterward.
        """
        if (not self.is_fset) or (not other.is_fset):
            raise ValueError("both pcs should be with is_fset=True to multiply")
        # transform self and other to fsets
        self_old_fset, other_old_fset = self.is_fset, other.is_fset
        self.is_fset = True
        other.is_fset = True
        # multiply
        new_pc_dict = {}
        for partition1, prob1 in self._pc_dict.items():
            for partition2, prob2 in other._pc_dict.items():
                prob_mul = self.prob_class.fast_prod(prob1, prob2)
                new_partition = cartesian_product_partition(partition1, partition2)
                if new_partition in new_pc_dict:
                    new_pc_dict[new_partition] += prob_mul
                else:
                    new_pc_dict[new_partition] = prob_mul
        # check for zero event
        if self.prob_class.need_to_check_zero_events:
            new_pc_dict = {
                partition: prob
                for partition, prob in new_pc_dict.items()
                if not self.prob_class.is_zero(prob)
            }
        # return the old fset option
        self.is_fset = self_old_fset
        other.is_fset = other_old_fset
        # return new pc
        return PartitionsCollection(
            pc_dict=new_pc_dict, prob_class=self.prob_class, is_fset=True
        )

    def items(self):
        """return the items of self.pc_dict"""
        return self._pc_dict.items()

    def add(self, other: Union[Dict, "PartitionsCollection"]):
        """
        Adds the partitions and probabilities from another collection or dictionary.

        Parameters
        ----------
        other : Union[Dict, PartitionsCollection]
            The partitions to add. Can be a dictionary or another PartitionsCollection.

        Raises
        ------
        ValueError
            If `other` is a PartitionsCollection and its `is_fset` attribute
            does not match the current collection.

        Notes
        -----
        - If `other` is a dictionary, its keys should match the format of this collection.
        - Modifies the current object in place.
        """
        if isinstance(other, PartitionsCollection) and (
            self._is_fset != other._is_fset
        ):
            raise ValueError("both collections should have the same is_fset value")
        # TODO: validate input
        add_dict_to_dict(self._pc_dict, other, inplace=True)

    def __add__(
        self, other: Union[Dict, "PartitionsCollection"]
    ) -> "PartitionsCollection":
        if isinstance(other, PartitionsCollection) and (
            self._is_fset != other._is_fset
        ):
            raise ValueError("both collections should have the same is_fset value")
        new_pc_dict = add_dict_to_dict(self._pc_dict, other._pc_dict, inplace=False)
        return PartitionsCollection(
            pc_dict=new_pc_dict, prob_class=self.prob_class, is_fset=self._is_fset
        )

    def __eq__(self, other: "PartitionsCollection") -> bool:
        if not isinstance(other, PartitionsCollection):
            return False
        if self.is_fset != other.is_fset:
            return False
        if len(self.pc_dict) != len(other.pc_dict):
            return False
        for partition, self_prob in self.items():
            if partition not in other.pc_dict:
                return False
            other_prob = other.pc_dict[partition]
            if np.abs(self_prob - other_prob) > 1e-7:
                return False
        return True

    def multiply_prob(self, prob) -> "PartitionsCollection":
        """
        Multiplies all probabilities in the collection by a given probability value.

        Parameters
        ----------
        prob : float
            The probability value to multiply with each partition's probability.

        Returns
        -------
        PartitionsCollection
            A new PartitionsCollection with updated probabilities.

        Notes
        -----
        - Removes partitions with zero probabilities after multiplication.
        - Uses `fast_prod` for `FSet` representations; otherwise, uses simple multiplication.
        """
        # create new_pc_dict
        new_pc_dict = {}
        for partition, _prob in self._pc_dict.items():
            multiply = self.prob_class.fast_prod if self.is_fset else lambda x, y: x * y
            prob_mul = multiply(_prob, prob)
            if not self.prob_class.is_zero(prob_mul):
                new_pc_dict[partition] = prob_mul
        return PartitionsCollection(
            pc_dict=new_pc_dict, prob_class=self.prob_class, is_fset=self._is_fset
        )

    def _raise_implemented_only_for_fset_pc(self, method_name):
        if not self.is_fset:
            raise ValueError(
                f"Method {method_name} is implemented only for pcs where pc.is_fset==True"
            )

    def minimal_form(self) -> "PartitionsCollection":
        self._raise_implemented_only_for_fset_pc("minimal_form")
        new_pc_dict = {}
        for partition, prob in new_pc_dict.items():
            new_partition = FSet.to_minimal_form(partition)
            add_dict_to_dict(new_pc_dict, {new_partition: prob})
        return PartitionsCollection(
            pc_dict=new_pc_dict, prob_class=self.prob_class, is_fset=self.is_fset
        )

    def zeta_trasform(self):
        if self._is_zeta:
            return self.copy()

    def create_lattice(self) -> nx.DiGraph:
        lattice = nx.DiGraph()
        for partition in self.pc_dict.keys():
            if (partition in lattice.nodes) and lattice.out_degree(partition) > 0:
                continue
            partition_list = [list(comp) for comp in partition]
            for i_comp, comp in enumerate(partition_list):
                for i_elm, elm in enumerate(comp):
                    new_partition = [
                        comp_p if i_comp != j_comp else
                        comp_p[:i_elm] + comp_p[i_elm + 1:]
                        for j_comp, comp_p in enumerate(partition)
                    ]
                    new_partition = remove_singleton(new_partition)
                    if new_partition in self.pc_dict:
                        lattice.add_edge(partition, new_partition)
        return lattice

MAX_LEN = 0
MAX_X = None


def cartesian_product_partitions_collection_list(
    partitions_collection_list: List[PartitionsCollection],
) -> PartitionsCollection:
    """
    Computes the Cartesian product of a list of PartitionsCollection objects.

    Parameters
    ----------
    partitions_collection_list : List[PartitionsCollection]
        A list of PartitionsCollection objects to multiply.

    Returns
    -------
    PartitionsCollection
        A new PartitionsCollection representing the Cartesian product.

    Raises
    ------
    ValueError
        If the input list is empty.

    Notes
    -----
    - Operates in decreasing order of collection sizes for efficiency.
    """
    global MAX_LEN, MAX_X
    if len(partitions_collection_list) == 0:
        raise ValueError("partitions_collection_list is empty")
    parti_collection_list_sorted = sorted(
        partitions_collection_list, key=len, reverse=True
    )
    max_len = max(len(x.pc_dict) for x in parti_collection_list_sorted)
    if max_len > MAX_LEN:
        MAX_LEN = max_len
        MAX_X = [len(x.pc_dict) for x in parti_collection_list_sorted]
    parti_collection_res = parti_collection_list_sorted[0]
    for pc in parti_collection_list_sorted[1:]:
        parti_collection_res *= pc
    return parti_collection_res
