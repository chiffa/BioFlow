"""
These methods are responsible for generation of pairs of nodes for which we will be calculating
and summing the flow.
"""
from itertools import combinations, repeat, product, cycle
from copy import copy
import random
from typing import Union, List, Tuple
import collections.abc
import numpy as np

from bioflow.utils.log_behavior import get_logger
from bioflow.utils.general_utils import _is_int

log = get_logger(__name__)


def reduce_and_deduplicate_sample(sample: Union[List[int], List[Tuple[int, float]]]) \
        -> List[Tuple[int, float]]:
    """
    Deduplicates the nodes found in the sample by adding weights of duplicated nodes. In case a
    list of node ids only is provided, transforms them into a weighted list with all weights set
    to 1.

    :param sample: sample to deduplicate and/or add weights to
    :return:
    """

    if _is_int(sample[0]):
        sample = [(node_id, 1) for node_id in sample]


    np_sample = np.array(sample).astype(np.float)
    u, c = np.unique(np_sample[:, 0], return_counts=True)
    dup = u[c > 1]

    if len(dup) > 0:

        for dup_id in dup:
            _mask = np_sample[:, 0] == dup_id
            sum_weight = np.sum(np_sample[_mask, 1])
            np_sample[_mask, 1] = 0
            np_sample[_mask[0], 1] = sum_weight

        sample = np_sample.tolist()
        sample = [(node_id, _value) for node_id, _value in sample if _value > 0]

    return sample


def evaluate_ops(prim_len: int, sec_len: int,
                 sparse_rounds: int = -1) -> float:
    """
    Evaluates the number of total node pair flow computations needed to calculate the complete
    flow in the sample according to the general_flow policy.

    :param prim_len: length of the primary set
    :param sec_len: length of the secondary set
    :param sparse_rounds: sparse rounds.
    :return:
    """

    if sparse_rounds < 1:

        if sec_len == 0:
            return prim_len**2 / 2

        else:
            return prim_len * sec_len

    else:

        if sec_len == 0:
            return prim_len * sparse_rounds // 2

        elif sec_len == 1:
            return prim_len

        else:
            return prim_len * sparse_rounds


def reduce_ops(prim_len, sec_len, max_ops) -> int:
    """
    Determines the sparse_rounds parameter that needs to be used in order to maintain the total
    number of pairs of nodes needed to calculate the complete flow in the sample according to the
    general_flow_policy.

    :param prim_len: length of the primary set
    :param sec_len: length of the secondary set
    :param max_ops: maximum allowed number of node pairs
    :return:
    """

    if evaluate_ops(prim_len, sec_len) < max_ops:
        return -1

    else:
        if sec_len == 0:
            return max(max_ops // (prim_len//2), 5)
        elif sec_len == 1:
            return -1
        else:
            return max(max_ops // prim_len, 5)


def general_flow(sample: Union[List[int], List[Tuple[int, float]]],
                 secondary_sample: Union[List[int], List[Tuple[int, float]], None] = None,
                 sparse_rounds: int = -1) -> List[Tuple[Tuple[int, float], Tuple[int, float]]]:
    """
    Performs the information flow computation best matching the provided parameters.

    :param sample: primary sample of nodes
    :param secondary_sample: secondary sample of nodes
    :param sparse_rounds: sparse rounds, in case samples are too big
    :return:
    """

    # TODO: what if we have an overlap between the items in the primary and the secondary
    #  samples?

    sample = reduce_and_deduplicate_sample(sample)

    if secondary_sample is None:  # connex flow

        if sparse_rounds > 0:
            list_of_pairs = []
            for _ in repeat(None, sparse_rounds):
                idx_list_c = copy(sample)
                random.shuffle(idx_list_c)
                list_of_pairs += list(zip(idx_list_c[:len(idx_list_c) // 2],
                                     idx_list_c[len(idx_list_c) // 2:]))
        else:
            list_of_pairs = [(i, j) for i, j in combinations(set(sample), 2)]

    else:

        secondary_sample = reduce_and_deduplicate_sample(secondary_sample)

        if sparse_rounds > 0:
            list_of_pairs = []
            for _ in range(sparse_rounds):
                idx_list_1 = copy(sample)
                idx_list_2 = copy(secondary_sample)
                random.shuffle(idx_list_1)
                random.shuffle(idx_list_2)
                list_of_pairs += list(zip(idx_list_1, cycle(idx_list_2)))

        else:
            list_of_pairs = [(i, j) for i, j in product(sample, secondary_sample)]

    return list_of_pairs

