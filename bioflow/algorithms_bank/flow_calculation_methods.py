"""
These methods are responsible for generation of pairs of nodes for which we will be calculating
and summing the flow.
"""
from itertools import combinations, repeat, product
from copy import copy
import random
from typing import Union, List, Tuple
import collections.abc
import numpy as np

from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)


def _is_int(_obj):
    try:
        int(_obj)
    except TypeError or ValueError as e:
        return False
    else:
        return True


def reduce_and_deduplicate_sample(sample: Union[List[int], List[Tuple[int, float]]]) \
        -> List[Tuple[int, float]]:

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
            np_sample[_mask[0], 1] = sum_weight  # INTEST: check that the _mask[0] does not fail

        sample = np_sample.tolist()
        sample = [(node_id, _value) for node_id, _value in sample if _value > 0]

    return sample


def general_flow(sample: Union[List[int], List[Tuple[int, float]]],
                 secondary_sample: Union[List[int], List[Tuple[int, float]], None] = None,
                 sparse_rounds: int = -1) -> List[Tuple[Tuple[int, float],
                                                        Tuple[int, float]]]:

    sample = reduce_and_deduplicate_sample(sample)

    if secondary_sample is None:  # connex flow

        if sparse_rounds > 0:
            list_of_pairs = []
            for _ in repeat(None, sparse_rounds):
                idx_list_c = copy(sample)  # CURRENTPASS: move out?
                random.shuffle(idx_list_c)
                list_of_pairs += list(zip(idx_list_c[:len(idx_list_c) // 2],
                                     idx_list_c[len(idx_list_c) // 2:]))
        else:
            list_of_pairs = [(i, j) for i, j in combinations(set(sample), 2)]

    else:

        secondary_sample = reduce_and_deduplicate_sample(secondary_sample)

        if sparse_rounds > 0:
            list_of_pairs = []
            for _ in repeat(None, sparse_rounds):  # CURRENTPASS range?
                idx_list_1 = copy(sample)  # CURRENTPASS: move out?
                idx_list_2 = copy(secondary_sample)  # CURRENTPASS: move out?
                random.shuffle(idx_list_1)
                random.shuffle(idx_list_2)
                list_of_pairs += list(zip(idx_list_1, idx_list_2))

        else:
            list_of_pairs = [(i, j) for i, j in product(sample, secondary_sample)]

    return list_of_pairs

