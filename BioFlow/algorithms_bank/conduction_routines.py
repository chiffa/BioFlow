"""
Module containing the the general routines for processing of conduction matrices with
IO current arrays.
"""
import random
from copy import copy
import numpy as np
from matplotlib import pyplot as plt
from itertools import combinations, repeat
from scipy.sparse import csc_matrix, diags, triu, lil_matrix
from scipy.sparse.linalg import eigsh
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky
from BioFlow.utils.log_behavior import logger as log
from BioFlow.internals_config import fudge
from BioFlow.utils.Linalg_routines import cluster_nodes, average_off_diag_in_sub_matrix, \
    average_interset_linkage, normalize_laplacian


def sparse_abs(sparse_matrix):
    """
    Recovers an absolute value of a sparse matrix

    :param sparse_matrix: sparse matrix for which we want to recover the absolute
    :return: absolute of that matrix
    """
    sign = sparse_matrix.sign()
    return sparse_matrix.multiply(sign)


def get_potentials_from_solver(laplacian_solver, io_array):
    """
    Solver wrapper for code clarity

    :param laplacian_solver: laplacian system solver
    :param io_array: array of currents for each node in system
    :return: potential in each node
    """
    return laplacian_solver(io_array)


def build_sink_source_current_array(io_index_pair, shape):
    """
    converts index pair to a solver-compatible array

    :param shape: shape of the conductance matrix
    :param io_index_pair: pair of indexes where sinks/source pair is
    """
    io_array = np.zeros((shape[0], 1))
    io_array[io_index_pair[0], 0], io_array[io_index_pair[1], 0] = (1.0, -1.0)
    return io_array


def get_potentials(conductivity_laplacian, io_index_pair):
    """
    Recovers voltages based on the conductivity Laplacian and the IO array

    :param conductivity_laplacian:
    :param io_index_pair:

    :return: array of potential in each node
    """
    # TODO: technically, Cholesky is not the best solver => change it if needed
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)
    io_array = build_sink_source_current_array(
        io_index_pair, conductivity_laplacian.shape)
    return get_potentials_from_solver(solver, io_array)


def get_current_matrix(conductivity_laplacian, voltages):
    """
    Recovers the current matrix based on the conductivity laplacian and voltages in each node.

    :param conductivity_laplacian:
    :param voltages:

    :return: matrix where M[i,j] = current intesity from i to j. Antisymmetric,
    if current i->j : coeff is positive, else coeff is negative.
    :return: matrix where M[i,j] = current intesity from i to j. Triangular superior,
    if current is from j to i, term is positive, otherwise it is negative.
    :rtype: scipy.sparse.lil_matrix
    """
    diag_voltages = lil_matrix(diags(voltages.T.tolist()[0], 0))
    corr_conductance_matrix = conductivity_laplacian - \
        lil_matrix(diags(conductivity_laplacian.diagonal(), 0))
    currents = diag_voltages.dot(
        corr_conductance_matrix) - corr_conductance_matrix.dot(diag_voltages)
    return currents, triu(currents)


def get_current_through_nodes(non_redundant_current_matrix):
    """

    :param non_redundant_current_matrix: non-redundant (i.e. triangular superior) matrix of
    currents through a conduction system
    :return : current through the individual nodes based on the current matrix as defined in
     the get_current_matrix module
    :rtype: numpy.array
    """
    poscurr = lil_matrix(non_redundant_current_matrix.shape)
    poscurr[
        non_redundant_current_matrix > 0.0] = non_redundant_current_matrix[
        non_redundant_current_matrix > 0.0]
    negcurr = lil_matrix(non_redundant_current_matrix.shape)
    negcurr[
        non_redundant_current_matrix < 0.0] = non_redundant_current_matrix[
        non_redundant_current_matrix < 0.0]
    s = np.array(poscurr.sum(axis=1).T - negcurr.sum(axis=0))
    r = np.array(poscurr.sum(axis=0) - negcurr.sum(axis=1).T)
    ret = copy(s)
    ret[r > s] = r[r > s]
    ret = list(ret.flatten())
    return ret


def get_pairwise_flow(
        conductivity_laplacian,
        index_list,
        cancellation=False,
        potential_dominated=True):
    """
    Performs a pairwaise computation and summation of the

    :param potential_dominated: if set to True, the computation is done by injecting constant
    potential difference into the system, not a constant current.
    :param conductivity_laplacian:  Laplacian representing the conductivity
    :param index_list: list of the indexes acting as current sources/sinks
    :param cancellation: if True, the conductance values that are induced by single nodes
    interactions will be cancelled
    :return: current matrix for the flow system; current through each node.
    """
    current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)

    for counter, (i, j) in enumerate(combinations(index_list, 2)):
        log.debug('getting pairwise flow %s out of %s', counter, len(index_list) ** 2 / 2)
        io_array = build_sink_source_current_array(
            (i, j), conductivity_laplacian.shape)
        voltages = get_potentials_from_solver(solver, io_array)
        currents_full, current_upper = get_current_matrix(
            conductivity_laplacian, voltages)

        if potential_dominated:
            potential_diff = abs(voltages[i, 0] - voltages[j, 0])
            current_upper = current_upper / potential_diff

        current_accumulator += sparse_abs(current_upper)

    if cancellation:
        sinks_no = len(index_list)
        current_accumulator /= (sinks_no * (sinks_no - 1) / 2)

    return current_accumulator


def get_better_pairwise_flow(
        conductivity_laplacian,
        index_list,
        cancellation=True,
        memory_source=None):
    """
    Performs a pairwise computation and summation of the pairwise_flow

    :param memory_source: dictionary of memoized tension and current flow through the circuit
    :param conductivity_laplacian:  Laplacian representing the conductivity
    :param index_list: list of the indexes acting as current sources/sinks
    :param cancellation: if True, the conductance values that are induced by single nodes
    interactions will be cancelled.
    :return: current matrix for the flow system
    :return: current through each node.
    """
    up_pair_2_voltage_current = {}
    current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)

    for counter, (i, j) in enumerate(combinations(set(index_list), 2)):
        log.debug('getting pairwise flow %s out of %s', counter, len(index_list) ** 2 / 2)

        if memory_source:
            potential_diff, current_upper = memory_source[
                tuple(sorted((i, j)))]

        else:
            io_array = build_sink_source_current_array(
                (i, j), conductivity_laplacian.shape)
            voltages = get_potentials_from_solver(solver, io_array)
            _, current_upper = get_current_matrix(
                conductivity_laplacian, voltages)
            potential_diff = abs(voltages[i, 0] - voltages[j, 0])
            up_pair_2_voltage_current[
                tuple(
                    sorted(
                        (i, j)))] = (
                potential_diff, current_upper)

        if potential_diff != 0:
            current_upper = current_upper / potential_diff
        else:
            print i, j
            raise Warning('Potential difference is nul')
        current_accumulator += sparse_abs(current_upper)

    if cancellation:
        ln = len(index_list)
        current_accumulator /= (ln * (ln - 1) / 2)

    return current_accumulator, up_pair_2_voltage_current


def sample_pairwise_flow(
        conductivity_laplacian,
        index_list,
        resamples,
        cancellation=False):
    """
    Performs sampling of pairwise flow in a conductance system.

    :param conductivity_laplacian: Laplacian representing the conductivity
    :param index_list: list of the indexes acting as current sources/sinks
    :param cancellation: if True, conductance values that are induced by single
    nodes interactions will be cancelled.
    :param resamples: number of times each element in idxlist will be sample.
    A reasonable minimal is such that len(idxlist)*resamples < 20 000

    :return: current matrix representing the flows from one node to the other. This
    flow is absolute and does not respect the Kirchoff's laws. However, it can be used to
    see the most important connections between the GO terms or Interactome and can be used to
    compute the flow through the individual nodes.
    """
    current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)
    list_of_pairs = []

    for _ in repeat(None, resamples):
        idx_list_c = copy(index_list)
        random.shuffle(idx_list_c)
        list_of_pairs += zip(idx_list_c[:len(idx_list_c) / 2], idx_list_c[len(idx_list_c) / 2:])

    for i, j in list_of_pairs:
        io_array = build_sink_source_current_array(
            (i, j), conductivity_laplacian.shape)
        voltages = get_potentials_from_solver(solver, io_array)
        currents_full, current_upper = get_current_matrix(
            conductivity_laplacian, voltages)
        current_accumulator += sparse_abs(current_upper)

    if cancellation:
        ln = len(index_list)
        current_accumulator /= (ln / 2 * resamples)

    return current_accumulator


def laplacian_reachable_filter(laplacian, reachable_indexes):
    """
    Modifies a laplacian matrix in a way so that only the elements that are reachable in the
    current iteration of flow computation are taken in account. Intended for the CO system
    computation.

    :warning: The only current alternative is usage of LU instead of cholesky, which is
    computationally more difficult and also requires reach-dependent computation to get
    an in and out flow to different GO terms

    :warning: An alternative is the construction of the individual laplacian
    for each new application

    :param laplacian: intial laplacian of directionless orientation
    :param reachable_indexes: indexes that are reachable from
    the nodes for which we want to perform the computation.
    :return: laplacian where all the lines and columns for terms that are not reachable are null.
    """
    arr = [0] * laplacian.shape[0]
    for index in reachable_indexes:
        arr[index] = 1
    diag_mat = diags(arr, 0, format="lil")
    re_laplacian = copy(laplacian)
    re_laplacian = diag_mat.dot(re_laplacian.dot(diag_mat))
    re_laplacian = re_laplacian - \
        diags(re_laplacian.diagonal(), 0, format="lil")
    d = (-re_laplacian.sum(axis=0)).tolist()[0]
    re_laplacian = re_laplacian + diags(d, 0, format="lil")
    return re_laplacian


def get_current_with_reach_limitations(
        inflated_laplacian, idx_pair, reach_limiter):
    """
    Recovers the current passing through a conduction system while enforcing the limitation
    on the directionality of induction of the GO terms

    :param inflated_laplacian: Laplacian containing the UP-GO relations in addition to
    purely GO-GO relations
    :param idx_pair: pair of indexes between which we want to compute the information flow
    :param reach_limiter: list of indexes to which we want to limit the reach
    :return:
    """
    reduced_laplacian = laplacian_reachable_filter(
        inflated_laplacian, reach_limiter)
    voltages = get_potentials(reduced_laplacian, (idx_pair[0], idx_pair[1]))
    _, current_upper = get_current_matrix(reduced_laplacian, voltages)
    potential_diff = abs(voltages[idx_pair[0], 0] - voltages[idx_pair[1], 0])
    current_upper = current_upper / potential_diff

    return current_upper, potential_diff


def perform_clustering(internode_tension, clusters, show=True):
    """
    Performs a clustering on the voltages of the nodes,

    :param internode_tension:
    :param clusters:
    :param show:
    """
    index_group = list(
        set([item for key in internode_tension.iterkeys() for item in key]))
    local_index = dict((UP, i) for i, UP in enumerate(index_group))
    rev_idx = dict((i, UP) for i, UP in enumerate(index_group))
    relmat = lil_matrix((len(index_group), len(index_group)))

    for (UP1, UP2), tension in internode_tension.iteritems():
        relmat[local_index[UP1], local_index[UP2]] = -1.0 / tension
        relmat[local_index[UP2], local_index[UP1]] = -1.0 / tension
        relmat[local_index[UP2], local_index[UP2]] += 1.0 / tension
        relmat[local_index[UP1], local_index[UP1]] += 1.0 / tension

    groups = cluster_nodes(relmat, clusters)

    relmat = normalize_laplacian(relmat)
    eigenvals, _ = eigsh(relmat)
    relmat = -relmat
    relmat.setdiag(1)

    groupsets = []
    group2average_offdiag = []
    for i in range(0, clusters):
        group_selector = groups == i
        group_idxs = group_selector.nonzero()[0].tolist()
        group2average_offdiag.append(
            (tuple(
                rev_idx[idx] for idx in group_idxs),
                len(group_idxs),
                average_off_diag_in_sub_matrix(
                relmat,
                group_idxs)))
        groupsets.append(group_idxs)

    remainder = average_interset_linkage(relmat, groupsets)

    clustidx = np.array([item for itemset in groupsets for item in itemset])
    relmat = relmat[:, clustidx]
    relmat = relmat[clustidx, :]

    mean_corr_array = np.array([[items, mean_corr]
                                for _, items, mean_corr in group2average_offdiag])

    if show:
        plt.imshow(relmat.toarray(), cmap='jet', interpolation="nearest")
        plt.colorbar()
        plt.show()

    return np.array(
        group2average_offdiag), remainder, mean_corr_array, eigenvals
