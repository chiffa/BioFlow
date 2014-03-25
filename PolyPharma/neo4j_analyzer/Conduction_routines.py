__author__ = 'ank'

import numpy as np
from scipy.sparse import lil_matrix


def get_voltages(concuctivity_laplacian, IOs):
    """
    Recovers voltages based on the conductivity laplacion and the IO array

    :param concuctivity_laplacian:
    :param IOs:

    :return: array of potential in each node
    """
    pass


def get_current_matrix(conductivity_laplacian, voltages):
    """
    Recovers the current matrix based on the conductivity laplacian and voltages in each node.

    :param conductivity_laplacian:
    :param voltages:

    :return: matrix where M[i,j] = current intesity from i to j. Triangular superior, if current is from j to i,
                term is positive, otherwise it is negative.
    :rtype: scipy.sparse.lil_matrix
    """
    pass


def get_current_through_nodes(current_matrix):
    """

    :param current_matrix:

    :return : current through the individual nodes based on the current matrix as defined in the get_current_matrix module
    :rtype: numpy.array
    """
    poscurr = lil_matrix(current_matrix.shape)
    poscurr[current_matrix > 0.0] = current_matrix[current_matrix > 0.0]
    negcurr = lil_matrix(current_matrix.shape)
    negcurr[current_matrix < 0.0] = current_matrix[current_matrix < 0.0]
    s = np.array(poscurr.sum(axis = 1).T - negcurr.sum(axis = 0))
    r = np.array(poscurr.sum(axis = 0) - negcurr.sum(axis = 1).T)
    ret = s
    ret[r > s] = r[r > s]
    ret = list(ret.flatten())
    return ret


def get_pairwise_flow(conductivity_laplacian, idxlist):
    """

    :param conductivity_laplacian:
    :param idxlist:

    :return: current matrix for the flow system
    :return: current through each node.
    """
    pass