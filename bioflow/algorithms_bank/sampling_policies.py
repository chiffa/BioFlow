"""
This module defines the policies that will be used in order to sample the information flow
patterns to compare with.

The general approach is a function that takes in any eventual parameters and outputs a list of
pairs of DB_Ids for which the flow will be calculated.
"""
import random
import hashlib
import json
import numpy as np
from typing import Union, List, Tuple
import collections.abc

from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def _is_int(_obj):
    try:
        int(_obj)
    except TypeError or ValueError as e:
        return False
    else:
        return True


def matched_sample_distribution(floats_arr: np.array, samples: int, granularity: int = 100):
    hist, bin_edges = np.histogram(floats_arr, bins=granularity, density=True)
    pad = np.arange(granularity)
    locations = np.choice(pad, samples, p=hist)

    samples = []

    for i in locations:
        samples.append(np.random.uniform(bin_edges[i], bin_edges[i+1]))

    return samples


def reduce_distribution(floats_arr: np.array):
    """
    Basically gets a distribution in the [0, 1] in 100 bins, rounds to the nearest 0.01

    :param floats_arr:
    :return:
    """
    normalized_arr = floats_arr / np.max(floats_arr)
    hist, bin_edges = np.histogram(normalized_arr, bins=100, density=True)
    rounded_hist = np.array(hist * 100).astype(np.int)

    return rounded_hist


def characterize_set(sample: Union[List[int], List[Tuple[int, float]]]):

    if sample is None:
        return 0, []

    if len(sample) == 1:
        return 1, []

    if _is_int(sample[0]):
        rounded_hist = [1]*100
        rounded_hist = np.array(rounded_hist).astype(np.int)
    else:
        rounded_hist = reduce_distribution(np.array(sample).astype(np.float)[:, 1])

    return len(sample), rounded_hist


def characterize_flow_parameters(sample: Union[List[int], List[Tuple[int, float]]],
                                 secondary_sample: Union[List[int], List[Tuple[int, float]], None],
                                 sparse_rounds: int):

    prim_len, prim_hist = characterize_set(sample)
    sec_len, sec_hist = characterize_set(secondary_sample)

    hash = hashlib.md5(json.dumps([prim_len, prim_hist,
                                   sec_len, sec_hist,
                                   sparse_rounds],
                                  sort_keys=True).encode('utf-8')).hexdigest()

    log.info('hashed a flow parameters from:\n'
             '%d/%s; \n'
             '%d/%s; \n'
             '%d \n'
             'to %s' % (prim_len, prim_hist, sec_len, sec_hist, sparse_rounds, hash))

    return prim_len, prim_hist, sec_len, sec_hist, sparse_rounds, hash


def matched_sampling(sample, secondary_sample, background, samples):

    if secondary_sample is None:
        pass

    pass


def set_signature(set_to_sign):
    if isinstance(set_to_sign[0], tuple):
        # it is an [(int, float)] weighted samples list, we hist sign floats
        floats = sorted([int(_item[1]*1000) for _item in set_to_sign])

    elif isinstance(set_to_sign[0], float):
         floats = sorted([int(_item*1000) for _item in set_to_sign])

    else:
        raise Exception('not sure how to sign set with first element', set_to_sign[0])

    return hashlib.md5(json.dumps(floats, sort_keys=True).encode('utf-8')).hexdigest()




def weighted_set_sampling(set_to_sample,
                          sample_size, samples,
                          weights=None):  # samples = iterations.
    if weights is None:
        for i in range(0, samples):
            random.shuffle(set_to_sample)
            yield i, set_to_sample[:sample_size]

    else:
        weights = np.array(weights)
        weights = weights / np.sum(weights)  # np. random choice expects probabilities summing to 1
        for i in range(0, samples):
            yield i, np.random.choice(set_to_sample, size=sample_size, replace=False, p=weights)

    pass


def double_weighted_set_sampling(set_to_sample, set_to_sample_2,
                                 sample_size, samples,
                                 weights1=None, weights_2=None):
    pass


def set_sampling(background, sample_size, iterations):
    """
    Basic sampling strategy for sets

    :param background:
    :param sample_size:
    :param iterations:
    :return:
    """
    for i in range(0, iterations):
        random.shuffle(background)
        yield i, background[:sample_size]


def degree_distribution_sampling():
    return []


def single_target_sampling():
    return []


def weight_distribution_sampling():
    return []


# CURRENTPASS: this needs to be replaced by a sampling_wrapper, that decide of the strategy based
#  on parameters provided.

active_default_sampling_policy = set_sampling