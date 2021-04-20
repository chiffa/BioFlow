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
    """
    Checks if an object is an int with a try-except loop

    :param _obj:
    :return:
    """
    try:
        int(_obj)
    except TypeError or ValueError as e:
        return False
    else:
        return True


def matched_sample_distribution(floats_arr: np.array, samples_no: int,
                                granularity: int = 100, logmode: bool = False) -> np.array:
    """
    Tries to guess a distribution of floats and sample from it.
    uses np.histogram with the number of bins equal to the granularity parameter. For each
    sample, selects which bin to sample and then picks from the bin a float according to a
    uniform distribution. if logmode is enabled, histogram will be in the log-space, as well as
    the sampling.

    :param floats_arr: array of floats for which to match the distribution
    :param samples_no: number of random samples to retrieve
    :param granularity: granularity at which to operate
    :param logmode: if sample in log-space
    :return: samples drawn from the empirically matched distribution
    """

    if logmode:
        floats_arr = np.log(floats_arr)   # INTEST: will crash if any are 0

    hist, bin_edges = np.histogram(floats_arr, bins=granularity, density=True)
    pad = np.arange(granularity)
    locations = np.choice(pad, samples_no, p=hist)

    samples = []

    for i in locations:
        samples.append(np.random.uniform(bin_edges[i], bin_edges[i+1]))

    if logmode:
        return np.exp(samples)

    else:
        return samples


def _reduce_distribution(floats_arr: np.array):
    """
    Basically gets a distribution in the [0, 1] in 100 bins, rounds to the nearest 0.01. Used for
    hashing and distribution matching

    :param floats_arr: floats for which to calculate the rounded distribution
    :return: rounded distribution
    """
    normalized_arr = floats_arr / np.max(floats_arr)
    bins = np.linspace(0, 1.001, 101)  # because floats round funny
    hist, bin_edges = np.histogram(normalized_arr, bins=bins, density=True)
    rounded_hist = np.array(hist * 100).astype(np.int)

    return rounded_hist


def _characterize_set(sample: Union[List[int], List[Tuple[int, float]]]):
    """
    None-robust helper function to characterize a sample set by its length, nature of items in
    teh sample and eventual distribution of weights within the sample.

    :param sample: sample to characterize
    :return: set length (0 if None), 1 if items are ids, 2 if ids and weights (0 if
    None), rounded distribution ([] if None or items are ids)
    """

    if sample is None:
        return 0, 0, []

    if len(sample) == 1:
        if _is_int(sample[0]):
            return 1, 1, []
        else:
            return 1, 2, []

    if _is_int(sample[0]):
        rounded_hist = [1] * 100
        rounded_hist = np.array(rounded_hist).astype(np.int)
        return len(sample), 1, rounded_hist.tolist()

    else:
        rounded_hist = _reduce_distribution(np.array(sample).astype(np.float)[:, 1])

        return len(sample), 2, rounded_hist.tolist()


def characterize_flow_parameters(sample: Union[List[int], List[Tuple[int, float]]],
                                 secondary_sample: Union[List[int], List[Tuple[int, float]], None],
                                 sparse_rounds: int):
    """
    Characterizes the primary and secondary sets and computes their hash, that can be used ot
    match similar samples for random sampling.

    :param sample: primary set
    :param secondary_sample: secondary set
    :param sparse_rounds: if sparse rounds are to be performed
    :return: first set length, shape, hist, second set length, shape, hist, sparse rounds, hash
    """

    prim_len, prim_shape, prim_hist = _characterize_set(sample)
    sec_len, sec_shape, sec_hist = _characterize_set(secondary_sample)

    _hash = hashlib.md5(json.dumps([prim_len, prim_shape, prim_hist,
                                   sec_len, sec_shape, sec_hist,
                                   sparse_rounds]).encode('utf-8')).hexdigest()

    log.debug('hashed a flow parameters from:\n'
             '%d/%d/%s; \n'
             '%d/%d/%s; \n'
             '%d \n'
             'to %s' % (prim_len, prim_shape, prim_hist,
                        sec_len, sec_shape, sec_hist,
                        sparse_rounds, _hash))

    return prim_len, prim_shape, prim_hist, sec_len, sec_shape, sec_hist, sparse_rounds, _hash


def _sample_floats(floats, float_sampling_method='exact', matched_distro_precision: int = 100):
    """
    A wrapper methods to sample a float distribution according to a method

    :param floats:
    :param float_sampling_method: exact (permutation of weights) | distro (trying to match the
        empirical distribution) | logdistro (trying to match the empirical distribution in the log
        space)
    :param matched_distro_precision: how closely to try to match the distribution (granularity
    parameter pass-through to the matched_sample_distribution)
    :return: sample of floats
    """

    if float_sampling_method == 'exact':
        ret_floats = floats.copy()
        np.random.shuffle(ret_floats)
        return ret_floats

    if float_sampling_method == 'distro':
        return matched_sample_distribution(floats, len(floats), granularity=matched_distro_precision)

    if float_sampling_method == 'logdistro':
        return matched_sample_distribution(floats, len(floats),
                                           granularity=matched_distro_precision, logmode=True)


def matched_sampling(sample, secondary_sample,
                     background, samples, float_sampling_method='exact'):
    """
    The general random sampling strategy that sample sets of the same size and shape as primary
    and secondary sample set and, if they are weighted, try to match the random sample weights
    according to the


    :param sample: primary sample set
    :param secondary_sample: secondary sample_set
    :param background: background of ids from which to sample
    :param samples: random samples wanted
    :param sampling_mode: exact/distro/logdistro. the sampling parametrization method ingesting
    all the parameters in a single string argument in the general case, here, a pass- through
    parameter for the _sample_floats function if samples are weighted and the distribution of
    weights is being matched.
    :return:
    """

    # CURRENTPASS: what if we have an overlap between the items in the primary and the secondary
    #  samples?

    if secondary_sample is None:

        if _is_int(sample[0]):  # INTEST: it will never be an int, but for safety ...
            for i in range(0, samples):
                random.shuffle(background)
                yield i, background[:len(sample)], None

        else:
            for i in range(0, samples):
                random.shuffle(background)
                id_loads = background[:len(sample)]
                float_part = _sample_floats(np.array(sample)[:, 1], float_sampling_method)
                ids_and_floats = [(_id, _float) for _id, _float in zip(id_loads, float_part)]
                yield i, ids_and_floats, None

    else:

        if _is_int(sample[0]):
            for i in range(0, samples):
                random.shuffle(background)
                yield i, background[:len(sample)], \
                      background[-len(secondary_sample):]

        else:
            for i in range(0, samples):
                random.shuffle(background)

                id_loads = background[:len(sample)]
                float_part = _sample_floats(np.array(sample)[:, 1], float_sampling_method)
                ids_and_floats = [(_id, _float) for _id, _float in zip(id_loads, float_part)]


                sec_id_loads = background[-len(secondary_sample):]
                sec_float_part = _sample_floats(np.array(secondary_sample)[:, 1], float_sampling_method)
                sec_ids_and_floats = [(_id, _float) for _id, _float
                                      in zip(sec_id_loads, sec_float_part)]

                yield i, ids_and_floats, sec_ids_and_floats
