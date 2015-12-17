"""
Module containing the the general routines for processing of conduction matrices with
IO current arrays.
"""
import random
from copy import copy
import numpy as np
from itertools import combinations, repeat
from scipy.sparse import csc_matrix, diags, triu, lil_matrix
from scipy.sparse.linalg import eigsh
# noinspection PyUnresolvedReferences
from scikits.sparse.cholmod import cholesky
import warnings
from BioFlow.utils.log_behavior import logger as log
from BioFlow.utils.dataviz import render_2d_matrix
from BioFlow.internal_configs import fudge
from BioFlow.utils.linalg_routines import cluster_nodes, average_off_diag_in_sub_matrix, \
    average_interset_linkage, normalize_laplacian

# TODO: we have to problems here: wrong solver and wrong laplacian
#   1) we are using a Cholesky solver on a system that by definition has at least one nul eigval
#   2) we are using the same laplacian matrix for all the calculations. However this is wrong:
#   we need to account for the fact that we are adding external sink/sources by adding 1
#   to the diagonal terms of the matrix that are being used as sinks/sources


def sparse_abs(sparse_matrix):
    """
    Recovers an absolute value of a sparse matrix

    :param sparse_matrix: sparse matrix for which we want to recover the absolute.
    :return: absolute of that matrix
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Changing the sparsity structure")

        sparse_matrix = csc_matrix(sparse_matrix)
        sign = sparse_matrix.sign()
        ret_mat = sparse_matrix.multiply(sign)

    return ret_mat


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
    return solver(io_array)


def get_current_matrix(conductivity_laplacian, node_potentials):
    """
    Recovers the current matrix based on the conductivity laplacian and voltages in each node.

    :param conductivity_laplacian:
    :param node_potentials:
    :return: matrix where M[i,j] = current intensity from i to j. Assymteric and Triangular
     superior iof the assymetric one. if current is from j to i, term is positive, otherwise
     it is negative.
    :rtype: scipy.sparse.lil_matrix
    """
    diag_voltages = lil_matrix(diags(node_potentials.T.tolist()[0], 0))
    corr_conductance_matrix = conductivity_laplacian - \
        lil_matrix(diags(conductivity_laplacian.diagonal(), 0))
    currents = diag_voltages.dot(corr_conductance_matrix) - \
        corr_conductance_matrix.dot(diag_voltages)
    return currents, triu(currents)


def get_current_through_nodes(triu_current_matrix):
    """
    Recovers current flowing through each node

    :param triu_current_matrix: non-redundant (i.e. triangular superior) matrix of
    currents through a conduction system
    :return : current through the individual nodes based on the current matrix as defined in
     the get_current_matrix module
    :rtype: numpy.array
    """
    pos_curr = lil_matrix(triu_current_matrix.shape)
    pos_curr[
        triu_current_matrix > 0.0] = \
        triu_current_matrix[triu_current_matrix > 0.0]
    neg_curr = lil_matrix(triu_current_matrix.shape)
    neg_curr[
        triu_current_matrix < 0.0] = \
        triu_current_matrix[triu_current_matrix < 0.0]
    s = np.array(pos_curr.sum(axis=1).T - neg_curr.sum(axis=0))
    r = np.array(pos_curr.sum(axis=0) - neg_curr.sum(axis=1).T)
    ret = copy(s)
    ret[r > s] = r[r > s]
    ret = list(ret.flatten())
    return ret


def laplacian_reachable_filter(laplacian, reachable_indexes):
    """
    Transforms a matrix to make sure only reachable elements are kept.

    The only current alternative is usage of LU instead of Cholesky, which is
    computationally more difficult and also requires reach-dependent computation to get
    an in and out flow to different GO terms

    An alternative is the construction of the individual laplacian
    for each new application

    :param laplacian: initial laplacian of directionless orientation
    :param reachable_indexes: indexes that are reachable from
    the nodes for which we want to perform the computation.
    :return: laplacian where all the lines and columns for terms that are not reachable are null.
    """
    pad_array = [0] * laplacian.shape[0]
    for index in reachable_indexes:
        pad_array[index] = 1
    diagonal_pad = diags(pad_array, 0, format="lil")
    re_laplacian = copy(laplacian)
    re_laplacian = diagonal_pad.dot(re_laplacian.dot(diagonal_pad))
    re_laplacian = re_laplacian - \
        diags(re_laplacian.diagonal(), 0, format="lil")
    d = (-re_laplacian.sum(axis=0)).tolist()[0]
    re_laplacian = re_laplacian + diags(d, 0, format="lil")
    return re_laplacian


def edge_current_iteration(conductivity_laplacian, index_pair,
                           solver=None, reach_limiter=None):
    """
    Master edge current retriever

    :param conductivity_laplacian:
    :param solver:
    :param index_pair:
    :param reach_limiter:
    :return: potential_difference, triu_current
    """
    if solver is None and reach_limiter is None:
        log.warning('edge current computation could be accelerated by using a shared solver')

    if reach_limiter:
        conductivity_laplacian = laplacian_reachable_filter(conductivity_laplacian, reach_limiter)
        solver = None

    i, j = index_pair
    io_array = build_sink_source_current_array((i, j), conductivity_laplacian.shape)
    if solver:
        voltages = solver(io_array)
    else:
        voltages = get_potentials(conductivity_laplacian, (i, j))
    _, current_upper = get_current_matrix(conductivity_laplacian, voltages)
    potential_diff = abs(voltages[i, 0] - voltages[j, 0])

    return potential_diff, current_upper


def master_edge_current(conductivity_laplacian, index_list,
                        cancellation=True, potential_dominated=True, sampling=False,
                        sampling_depth=10, memory_source=None, memoization=None):
    """
    master method for all the required edge current calculations

    :param conductivity_laplacian:
    :param index_list:
    :param cancellation:
    :param potential_dominated:
    :param sampling:
    :param sampling_depth:
    :param memory_source:
    :param memoization:
    :return:
    """
    # TODO: remove memoization in order to reduce overhead nad mongo database load

    # generate index list in agreement with the sampling strategy
    if sampling:
        list_of_pairs = []
        for _ in repeat(None, sampling_depth):
            idx_list_c = copy(index_list)
            random.shuffle(idx_list_c)
            list_of_pairs += zip(idx_list_c[:len(idx_list_c) / 2],
                                 idx_list_c[len(idx_list_c) / 2:])
    else:
        list_of_pairs = [(i, j) for i, j in combinations(set(index_list), 2)]

    total_pairs = len(list_of_pairs)

    up_pair_2_voltage_current = {}
    current_accumulator = lil_matrix(conductivity_laplacian.shape)
    solver = cholesky(csc_matrix(conductivity_laplacian), fudge)

    # run the main loop on the list of indexes in agreement with the memoization strategy:
    for counter, (i, j) in enumerate(list_of_pairs):
        log.debug('getting pairwise flow %s out of %s', counter + 1, total_pairs)

        if memory_source and tuple(sorted((i, j))) in memory_source.keys():
            potential_diff, current_upper = memory_source[tuple(sorted((i, j)))]

        else:
            potential_diff, current_upper = \
                edge_current_iteration(conductivity_laplacian, (i, j),
                                       solver=solver)

        if memoization:
            up_pair_2_voltage_current[tuple(sorted((i, j)))] = \
                (potential_diff, current_upper)

        # normalize to potential, if needed
        if potential_dominated:

            if potential_diff != 0:
                current_upper = current_upper / potential_diff

            # warn if potential difference is null or close to it
            else:
                log.warning('pairwise flow. On indexes %s %s potential difference is null. %s',
                            i, j, 'Tension-normalization was aborted')

        current_accumulator += sparse_abs(current_upper)

    if cancellation:
        current_accumulator /= total_pairs

    return current_accumulator, up_pair_2_voltage_current


def group_edge_current(conductivity_laplacian, index_list,
                       cancellation=False, potential_dominated=True):
    """
    Performs a pairwise computation and summation of the

    :param conductivity_laplacian:  Laplacian representing the conductivity
    :param index_list: list of the indexes acting as current sources/sinks
    :param cancellation: if True, conductance would be normalized to number of sinks used
    :param potential_dominated: if set to True, the computation is done by injecting constant
    potential difference into the system, not a constant current.
    :return: current matrix for the flow system; current through each node.
    """
    current_accumulator, _ = master_edge_current(
        conductivity_laplacian, index_list,
        cancellation=cancellation,
        potential_dominated=potential_dominated)

    return current_accumulator


def group_edge_current_memoized(conductivity_laplacian, index_list,
                                cancellation=True, memory_source=None):
    """
    Performs a pairwise computation and summation of the pairwise_flow

    :param conductivity_laplacian: Laplacian representing the conductivity
    :param index_list: list of the indexes acting as current sources/sinks
    :param cancellation: if True, conductance would be normalized to number of sinks used
    :param memory_source: dictionary of memoized tension and current flow through the circuit
    :return: current matrix for the flow system, current through each node.
    """
    return master_edge_current(conductivity_laplacian, index_list,
                               cancellation=cancellation,
                               memory_source=memory_source,
                               memoization=True)


def sample_group_edge_current(conductivity_laplacian, index_list, re_samples,
                              cancellation=False):
    """
    Performs sampling of pairwise flow in a conductance system.

    :param conductivity_laplacian: Laplacian representing the conductivity
    :param index_list: list of the indexes acting as current sources/sinks
    :param cancellation: if True, conductance would be normalized to number of sinks used
    :param re_samples: number of times each element in idxlist will be sample.
    A reasonable minimal is such that len(idxlist)*resamples < 20 000
    :return: current matrix representing the flows from one node to the other. This
    flow is absolute and does not respect the Kirchoff's laws. However, it can be used to
    see the most important connections between the GO terms or Interactome and can be used to
    compute the flow through the individual nodes.
    """

    current_accumulator, _ = master_edge_current(conductivity_laplacian, index_list,
                                                 cancellation=cancellation,
                                                 sampling=True,
                                                 sampling_depth=re_samples)

    return current_accumulator


def group_edge_current_with_limitations(inflated_laplacian, idx_pair, reach_limiter):
    """
    Recovers the current passing through a conduction system while enforcing the limitation
    on the directionality of induction of the GO terms

    :param inflated_laplacian: Laplacian containing the UP-GO relations in addition to
    purely GO-GO relations
    :param idx_pair: pair of indexes between which we want to compute the information flow
    :param reach_limiter: list of indexes to which we want to limit the reach
    :return:
    """
    inverter = edge_current_iteration(inflated_laplacian, idx_pair,
                                      reach_limiter=reach_limiter)

    return inverter[1]/inverter[0], inverter[0]


def perform_clustering(inter_node_tension, cluster_number, show=True):
    """
    Performs a clustering on the voltages of the nodes,

    :param inter_node_tension:
    :param cluster_number:
    :param show:
    """
    index_group = list(set([item
                            for key in inter_node_tension.iterkeys()
                            for item in key]))
    local_index = dict((UP, i) for i, UP in enumerate(index_group))
    rev_idx = dict((i, UP) for i, UP in enumerate(index_group))
    relations_matrix = lil_matrix((len(index_group), len(index_group)))

    for (UP1, UP2), tension in inter_node_tension.iteritems():
        # TODO: change the metric used to cluster the nodes.
        relations_matrix[local_index[UP1], local_index[UP2]] = -1.0 / tension
        relations_matrix[local_index[UP2], local_index[UP1]] = -1.0 / tension
        relations_matrix[local_index[UP2], local_index[UP2]] += 1.0 / tension
        relations_matrix[local_index[UP1], local_index[UP1]] += 1.0 / tension

    # underlying method is spectral clustering: do we really lie in a good zone for that?
    groups = cluster_nodes(relations_matrix, cluster_number)

    relations_matrix = normalize_laplacian(relations_matrix)
    eigenvals, _ = eigsh(relations_matrix)
    relations_matrix = -relations_matrix
    relations_matrix.setdiag(1)

    group_sets = []
    group_2_mean_off_diag = []
    for i in range(0, cluster_number):
        group_selector = groups == i
        group_indexes = group_selector.nonzero()[0].tolist()
        group_2_mean_off_diag.append(
            (tuple(
                rev_idx[idx] for idx in group_indexes),
                len(group_indexes),
                average_off_diag_in_sub_matrix(
                relations_matrix,
                group_indexes)))
        group_sets.append(group_indexes)

    remainder = average_interset_linkage(relations_matrix, group_sets)

    clustidx = np.array([item for itemset in group_sets for item in itemset])
    relations_matrix = relations_matrix[:, clustidx]
    relations_matrix = relations_matrix[clustidx, :]

    mean_corr_array = np.array([[items, mean_corr]
                                for _, items, mean_corr in group_2_mean_off_diag])

    if show:
        render_2d_matrix(relations_matrix.toarray(), 'Relationships matrix')

    return np.array(group_2_mean_off_diag), \
        remainder, \
        mean_corr_array, \
        eigenvals
