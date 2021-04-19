import numpy as np
from scipy.stats import gumbel_r
from typing import List

from bioflow.utils.log_behavior import get_logger


log = get_logger(__name__)


def get_neighboring_degrees(degree: int,
                            max_array: np.array,
                            nearest_degrees: int = 0,  # currently inactive
                            min_nodes: int = 10) -> List[float]:
    """
    Recovers the maximum flow achieved by nodes of a given degree for each run. On case the user
    requests it with nearest_degrees or min_nodes parameters, also recovers maximum flow
    values for the nodes of similar degrees or looks fro flow values in nearest degrees until
    at least `min_nodes` are found
    `nearest_degrees`

    :param degree: degree of the nodes
    :param max_array: maximum nodes for a given degree in each run
    :param nearest_degrees: the minimum number of the nearest gedgrees to look for
    :param min_nodes: the minimum number of nodes until which to look for neighbours
    :return:
    """

    max_set = max_array[:, max_array[1, :] == degree]
    max_set_red = max_set[0, :].tolist()

    if len(max_set_red) < min_nodes:
        temp_deg_plus = degree
        temp_deg_minus = degree

        while len(max_set_red) < min_nodes:
            temp_deg_minus -= 1
            temp_deg_plus += 1

            _max_set = max_array[:, max_array[1, :] == temp_deg_minus]
            max_set_red += _max_set[0, :].tolist()

            _max_set = max_array[:, max_array[1, :] == temp_deg_plus]
            max_set_red += _max_set[0, :].tolist()

    log.debug('deg: %d, list: %s' % (degree, max_set_red))

    return max_set_red


def get_p_val_by_gumbel(entry: np.array,
                        max_set_red: List[float]) -> np.array:
    """
    Recovers the statistical significance (p-value equivalent) by performing a gumbel test

    :param entry: the values achieved in the real hits information flow computation
    :param max_set_red: background set of maximum values achieved during blanc sparse_sampling runs
    :return:
    """

    params = gumbel_r.fit(max_set_red)
    mu = params[-2]
    beta = params[-1]

    frozen_gumbel = gumbel_r(loc=mu, scale=beta)
    a_95_low, a_95_high = frozen_gumbel.interval(0.95)

    log.debug('gumbel_r fit: mu %.2f, beta: %.2f, 95 alpha: %.2f, .%.2f'
              % (mu, beta, a_95_low, a_95_high))

    p_vals = 1 - frozen_gumbel.cdf(entry[0, :])

    return p_vals
