"""
These methods are responsible for generation of pairs of nodes for which we will be calculating
and summing the flow.
"""
from itertools import combinations, repeat
from copy import copy
import random
from typing import Union, List, Tuple


def connex_flow(sample: Union[List[int], List[Tuple[int, float]]],
                sparse: bool,
                sparse_rounds: int) -> Union[List[Tuple[int, int]],
                                             List[Tuple[Tuple[int, float],
                                                        Tuple[int, float]]]]:
    """
    Connex flow sampling routine

    :param sample:
    :param sparse:
    :param sparse_rounds:
    :return:
    """
    if sparse:
        list_of_pairs = []
        for _ in repeat(None, sparse_rounds):
            idx_list_c = copy(sample)
            random.shuffle(idx_list_c)
            list_of_pairs += list(zip(idx_list_c[:len(idx_list_c) // 2],
                                 idx_list_c[len(idx_list_c) // 2:]))
    else:
        list_of_pairs = [(i, j) for i, j in combinations(set(sample), 2)]

    return list_of_pairs


def star_flow(sample: Union[List[int], List[Tuple[int, float]]],
              center: Union[int, Tuple[int, float]],
              sparse: bool,
              spare_rounds: int) -> Union[List[Tuple[int, int]],
                                          List[Tuple[Tuple[int, float],
                                                     Tuple[int, float]]]]:
    return []


def biparty_flow(sample: Union[List[int], List[Tuple[int, float]]],
                 target_sample: Union[List[int], List[Tuple[int, float]]],
                 sparse: bool,
                 sparse_rounds: int) -> Union[List[Tuple[int, int]],
                                              List[Tuple[Tuple[int, float],
                                                         Tuple[int, float]]]]:
    return []


# TRACING: this needs to be injected into the main_flow_calc_loop
def flow_wrapper(sample: Union[List[int], List[Tuple[int, float]]],
                 center: Union[int, Tuple[int, float]] = None,
                 target_sample: Union[List[int], List[Tuple[int, float]]] = None,
                 sparse: bool = False,
                 sparse_rounds: int = -1) -> Tuple[str, Union[List[Tuple[int, int]],
                                                              List[Tuple[Tuple[int, float],
                                                                         Tuple[int, float]]]]]:

    # we return the name as well because that's what we use for storage as well
    # CURRENTPASS: not really, because we can store parameters that will decide the routine
    if center is not None:
        return 'star', star_flow(sample, center, sparse, sparse_rounds)

    if target_sample is not None:
        return 'biparty', biparty_flow(sample, target_sample, sparse, sparse_rounds)

    return 'connex', connex_flow(sample, sparse, sparse_rounds)

