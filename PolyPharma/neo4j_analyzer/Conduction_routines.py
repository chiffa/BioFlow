__author__ = 'ank'

import numpy as np
from scipy.sparse import lil_matrix, csc_matrix, diags, triu
from PolyPharma.configs import fudge
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky

def get_voltages_with_solver(laplacian_solver, IO_array):
    """
    A convinient solver wrapper to avoind calling all the time the Solver(IO_array) assignement routine

    :param laplacian_solver: object being able to solve the condctance laplacian system for the currents outside the system
    :param IO_array: array representing the currents between the nodes within the system and outside of the system
    :return: potential in each node
    """
    return laplacian_solver(IO_array)


def get_IO_currents_array(IO_index_pair,shape):
    """
    from a set of two indexes, builds and returns the np array that will be taken in by a solver to recover the
    voltages

    :param IO_index_pair: pair of indexes that represent the input/output nodes(indifferently)
    """
    IO_array = np.zeros((shape[0], 1))
    IO_array[IO_index_pair[0], 0], IO_array[IO_index_pair[1], 0] = (1.0, -1.0)
    return  IO_array


def get_voltages(conductivity_laplacian, IO_index_pair):
    """
    Recovers voltages based on the conductivity laplacion and the IO array

    :param concuctivity_laplacian:
    :param IOs:

    :return: array of potential in each node
    :return: difference of potentials between the source and the sink
    """
    Solver = cholesky(csc_matrix(conductivity_laplacian), fudge)
    IO_array = get_IO_currents_array(IO_index_pair)
    return get_voltages_with_solver(Solver, IO_array)


def get_current_matrix(conductivity_laplacian, voltages):
    """
    Recovers the current matrix based on the conductivity laplacian and voltages in each node.

    :param conductivity_laplacian:
    :param voltages:

    :return: matrix where M[i,j] = current intesity from i to j. Antisymmetric, if current i->j : coeff is positive, else
                coeff is negative.
    :return: matrix where M[i,j] = current intesity from i to j. Triangular superior, if current is from j to i,
                term is positive, otherwise it is negative.
    :rtype: scipy.sparse.lil_matrix
    """
    diag_Voltages = lil_matrix(diags(voltages.T.tolist()[0], 0))
    Corr_Conductance_Matrix = conductivity_laplacian - lil_matrix(diags(conductivity_laplacian.diagonal(), 0))
    currents = diag_Voltages.dot(Corr_Conductance_Matrix) - Corr_Conductance_Matrix.dot(diag_Voltages)
    return currents, triu(currents)


def get_current_through_nodes(non_redundant_current_matrix):
    """

    :param non_redundant_current_matrix: non-redundant (i.e. triangular superior) matrix of currents through a conduction
                                            system

    :return : current through the individual nodes based on the current matrix as defined in the get_current_matrix module
    :rtype: numpy.array
    """
    poscurr = lil_matrix(non_redundant_current_matrix.shape)
    poscurr[non_redundant_current_matrix > 0.0] = non_redundant_current_matrix[non_redundant_current_matrix > 0.0]
    negcurr = lil_matrix(non_redundant_current_matrix.shape)
    negcurr[non_redundant_current_matrix < 0.0] = non_redundant_current_matrix[non_redundant_current_matrix < 0.0]
    s = np.array(poscurr.sum(axis = 1).T - negcurr.sum(axis = 0))
    r = np.array(poscurr.sum(axis = 0) - negcurr.sum(axis = 1).T)
    ret = s
    ret[r > s] = r[r > s]
    ret = list(ret.flatten())
    return ret


def get_pairwise_flow(conductivity_laplacian, idxlist, cancellation=True):
    """

    :param conductivity_laplacian:
    :param idxlist:
    :param cancellation: if True, the conductance values that are induced by single nodes interactions will be cancelled.

    :return: current matrix for the flow system
    :return: current through each node.
    """
    pass


def sample_pairwise_flow(conductivity_laplacian, idxlist, resamples, cancellation=True):
    """
    Performs sampling of pairwise flow in a conductance system.

    :param conductivity_laplacian: Laplacian representing the conductivity
    :param idxlist: list of the indexes acting as current sources/sinks
    :param cancellation: if True, conductance values that are induced by single nodes interactions will be cancelled.
    :param resamples: number of times each element in idxlist will be sample. A reasonable minimal is such that len(idxlist)*resamples < 20 000

    :return: current matrix representing the flows from one node to the other. This flow is absolute and does not respect
                the Kirchoff's laws. However, it can be used to see the most important connections between the GO terms
                or Interactome and can be used to compute the flow through the individual nodes.
    """
    return 1