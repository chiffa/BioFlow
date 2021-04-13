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