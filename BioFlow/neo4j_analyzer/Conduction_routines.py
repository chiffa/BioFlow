"""
Module containing the the general routines for processing of conduction matrices with IO current arrays.
"""
__author__ = 'ank'

import random
import numpy as np
from scipy.sparse import csc_matrix, diags, triu
from BioFlow.configs2 import fudge
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky
from itertools import combinations, repeat
from copy import copy
from BioFlow.Utils.Linalg_routines import cluster_nodes, submatrix, remaineder_matrix, normalize_laplacian
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import eigsh
from matplotlib import pyplot as plt


def sparse_abs(sparse_matrix):
    """
    Recovers an absolute value of a sparse matrix

    :param sparse_matrix: sparse matrix for which we want to recover the absolute
    :return: absolute of that matrix
    """
    sign = sparse_matrix.sign()
    return sparse_matrix.multiply(sign)


def get_voltages_with_solver(laplacian_solver, IO_array):
    """
    A convinient solver wrapper to avoind calling all the time the Solver(IO_array) assignement routine

    :param laplacian_solver: object being able to solve the condctance laplacian system for the currents outside the system
    :param IO_array: array representing the currents between the nodes within the system and outside of the system
    :return: potential in each node
    """
    return laplacian_solver(IO_array)


def build_IO_currents_array(IO_index_pair,shape):
    """
    from a set of two indexes, builds and returns the np array that will be taken in by a solver to recover the
    voltages

    :param shape: shape of the conductance matrix
    :param IO_index_pair: pair of indexes that represent the input/output nodes(indifferently)
    """
    IO_array = np.zeros((shape[0], 1))
    IO_array[IO_index_pair[0], 0], IO_array[IO_index_pair[1], 0] = (1.0, -1.0)
    return  IO_array


def get_voltages(conductivity_laplacian, IO_index_pair):
    """
    Recovers voltages based on the conductivity laplacion and the IO array

    :param conductivity_laplacian:
    :param IO_index_pair:

    :return: array of potential in each node
    :return: difference of potentials between the source and the sink
    """
    Solver = cholesky(csc_matrix(conductivity_laplacian), fudge)
    IO_array = build_IO_currents_array(IO_index_pair,conductivity_laplacian.shape)
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
    ret = copy(s)
    ret[r > s] = r[r > s]
    ret = list(ret.flatten())
    return ret


def get_pairwise_flow(conductivity_laplacian, idxlist, cancellation=False, potential_dominated=True):
    """
    Performs a pairwaise computation and summation of the

    :param potential_dominated: if set to True, the computation is done by injecting constant potential difference into the
                            system, not a constant current.
    :param conductivity_laplacian:  Laplacian representing the conductivity
    :param idxlist: list of the indexes acting as current sources/sinks
    :param cancellation: if True, the conductance values that are induced by single nodes interactions will be cancelled.

    :return: current matrix for the flow system
    :return: current through each node.
    """
    Current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)

    for counter, (i, j) in enumerate(combinations(idxlist, 2)):
        # print 'getting pairwise flow %s out of %s' % (counter, len(idxlist)**2/2)
        IO_array = build_IO_currents_array((i,j), conductivity_laplacian.shape)
        voltages = get_voltages_with_solver(solver, IO_array)
        currents_full, current_upper = get_current_matrix(conductivity_laplacian, voltages)

        if potential_dominated:
            potential_diff = abs(voltages[i,0]-voltages[j,0])
            current_upper = current_upper/potential_diff

        Current_accumulator += sparse_abs(current_upper)

    if cancellation:
        ln = len(idxlist)
        Current_accumulator = Current_accumulator/(ln*(ln-1)/2)

    return Current_accumulator


def get_better_pairwise_flow(conductivity_laplacian, idxlist, cancellation=True, memoized=False, memory_source=None):
    """
    Performs a pairwaise computation and summation of the


    :param memory_source: dictionary of memoized tension and current flow through the circuit
    :param conductivity_laplacian:  Laplacian representing the conductivity
    :param idxlist: list of the indexes acting as current sources/sinks
    :param cancellation: if True, the conductance values that are induced by single nodes interactions will be cancelled.

    :return: current matrix for the flow system
    :return: current through each node.
    """
    UP_pair2voltage_current = {}
    Current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)

    for counter, (i, j) in enumerate(combinations(set(idxlist), 2)):

        # print 'getting pairwise flow %s out of %s' % (counter, len(idxlist)*(len(idxlist)-1)/2)

        if memory_source:
            potential_diff, current_upper = memory_source[tuple(sorted((i, j)))]

        else:
            IO_array = build_IO_currents_array((i, j), conductivity_laplacian.shape)
            voltages = get_voltages_with_solver(solver, IO_array)
            _, current_upper = get_current_matrix(conductivity_laplacian, voltages)
            potential_diff = abs(voltages[i, 0]-voltages[j, 0])
            UP_pair2voltage_current[tuple(sorted((i, j)))] = (potential_diff, current_upper)

        if potential_diff!=0:
            current_upper = current_upper/potential_diff
        else:
            print i, j
            raise Warning('Potential difference is nul')
        Current_accumulator += sparse_abs(current_upper)

    if cancellation:
        ln = len(idxlist)
        Current_accumulator = Current_accumulator/(ln*(ln-1)/2)

    return Current_accumulator, UP_pair2voltage_current



def sample_pairwise_flow(conductivity_laplacian, idxlist, resamples, cancellation=False):
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
    Current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)
    List_of_pairs = []

    for _ in repeat(None, resamples):
        L = copy(idxlist)
        random.shuffle(L)
        List_of_pairs += zip(L[:len(L)/2], L[len(L)/2:])

    for i, j in List_of_pairs:
        IO_array = build_IO_currents_array((i,j), conductivity_laplacian.shape)
        voltages = get_voltages_with_solver(solver, IO_array)
        currents_full, current_upper = get_current_matrix(conductivity_laplacian,voltages)
        Current_accumulator += sparse_abs(current_upper)

    if cancellation:
        ln = len(idxlist)
        Current_accumulator = Current_accumulator/(ln/2*resamples)

    return Current_accumulator


def laplacian_reachable_filter(laplacian, reachable_indexes):
    """
    Modifies a laplacian matrix in a way so that only the elemetns that are reachable in the current iteration of flow
    computation are taken in account. Intended for the CO system computation.

    :warning: The only current alternative is usage of LU instead of cholesky, which is computationally more difficult and also requires reach-dependent computation to
    get an in and out flow to different GO terms

    :warning: An alternative is the construction of the individual laplacian for each new application

    :param laplacian: intial laplacian of directionless orientation
    :param reachable_indexes: indexes that are reachable from the nodes for which we want to perform the computation.
    :return: laplacian where all the lines and colums for terms that are not reachable are null.
    """
    arr = [0]*laplacian.shape[0]
    for index in reachable_indexes:
        arr[index] = 1
    diag_mat = diags(arr, 0, format="lil")
    re_laplacian = copy(laplacian)
    re_laplacian = diag_mat.dot(re_laplacian.dot(diag_mat))
    re_laplacian = re_laplacian - diags(re_laplacian.diagonal(), 0, format="lil")
    d = (-re_laplacian.sum(axis = 0)).tolist()[0]
    re_laplacian = re_laplacian + diags( d, 0, format="lil")
    return re_laplacian



def get_current_with_reach_limitations(inflated_laplacian, Idx_pair, reach_limiter):
    """
    Recovers the current passing through a conduction system while enforcing the limitation on the directionality of induction of the GO terms

    :param inflated_laplacian: Laplacian containing the UP-GO relations in addition to purelyu GO-GO relations
    :param Idx_pair: pair of indexes between which we want to compute the information flow
    :param reach_limiter: list of indexes to which we want to limit the reach
    :return:
    """
    reduced_laplacian = laplacian_reachable_filter(inflated_laplacian, reach_limiter)
    voltages = get_voltages(reduced_laplacian, (Idx_pair[0], Idx_pair[1]))
    _, current_upper = get_current_matrix(reduced_laplacian, voltages)
    potential_diff = abs(voltages[Idx_pair[0], 0] - voltages[Idx_pair[1], 0])
    current_upper = current_upper / potential_diff

    return current_upper, potential_diff


def perform_clustering(internode_tension, clusters, show=True):
    """
    Performs a clustering on the voltages of the nodes,

    :param internode_tension:
    """
    idxgroup = list(set([ item for key in internode_tension.iterkeys() for item in key ]))
    local_index = dict((UP, i) for i, UP in enumerate(idxgroup))
    rev_idx = dict((i, UP) for i, UP in enumerate(idxgroup))
    relmat = lil_matrix((len(idxgroup), len(idxgroup)))

    for (UP1, UP2), tension in internode_tension.iteritems():
        relmat[local_index[UP1], local_index[UP2]] = -1.0/tension
        relmat[local_index[UP2], local_index[UP1]] = -1.0/tension
        relmat[local_index[UP2], local_index[UP2]] += 1.0/tension
        relmat[local_index[UP1], local_index[UP1]] += 1.0/tension

    groups = cluster_nodes(relmat, clusters)

    relmat = normalize_laplacian(relmat)
    eigenvals, _ = eigsh(relmat)
    relmat = -relmat
    relmat.setdiag(1)

    groupsets = []
    group2average_offdiag = []
    for i in range(0, clusters):
        group_selector = groups==i
        group_idxs = group_selector.nonzero()[0].tolist()
        group2average_offdiag.append((tuple(rev_idx[idx] for idx in group_idxs), len(group_idxs), submatrix(relmat, group_idxs)))
        groupsets.append(group_idxs)

    remainder = remaineder_matrix(relmat, groupsets)

    clustidx = np.array([item for itemset in groupsets for item in itemset])
    relmat = relmat[:, clustidx]
    relmat = relmat[clustidx, :]

    mean_corr_array = np.array([[items, mean_corr] for _, items, mean_corr in group2average_offdiag])

    if show:
        plt.imshow(relmat.toarray(), cmap='jet', interpolation="nearest")
        plt.colorbar()
        plt.show()

    return np.array(group2average_offdiag), remainder, mean_corr_array, eigenvals
