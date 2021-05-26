"""
Module containing the the general routines for processing of conduction matrices with
IO current arrays.
"""
import random
from copy import copy
from time import time
import datetime
import numpy as np
import importlib
from itertools import combinations, repeat
import scipy.sparse as spmat
# from scipy.sparse.linalg import eigsh
import scikits.sparse.cholmod as chmd
# from scikits.sparse.cholmod import cholesky, Factor
from scipy.sparse.linalg import splu
import warnings
from typing import Union, Tuple, List

from bioflow.utils.log_behavior import get_logger
# from bioflow.internal_configs import line_loss
from bioflow.algorithms_bank.flow_calculation_methods import general_flow
from bioflow.configs.main_configs import switch_to_splu, share_solver, memory_source_allowed, \
    node_current_in_debug, line_loss

log = get_logger(__name__)

# switch_to_splu = False
# # Looks like we are failing the normalization due to the matrix symmetry when using SPLU.
# # Which is expected - since we did simplifying assumptions about the Laplacian to be able to share it
#
# # In theory, we have to problems here: wrong solver and wrong laplacian
# #   1) we are using a Cholesky solver on a system that by definition has at least one nul eigval
# #   2) we are using the same laplacian matrix for all the calculations. However this is wrong:
# #   we need to account for the fact that we are adding external sink/sources by adding 1
# #   to the diagonal terms of the matrix that are being used as sinks/sources
# #
# #   A comparison with a proper SPLU showed we make an error that is extremely small but perform the calculation about 20x faster


def delete_row_csr(mat, i):
    """
    Deletes a row in a csr matrix (warning: inefficient)

    :param mat: matrix
    :param i: row index
    :return:
    """
    if not isinstance(mat, spmat.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    n = mat.indptr[i+1] - mat.indptr[i]
    if n > 0:
        mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
        mat.data = mat.data[:-n]
        mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
        mat.indices = mat.indices[:-n]
    mat.indptr[i:-1] = mat.indptr[i+1:]
    mat.indptr[i:] -= n
    mat.indptr = mat.indptr[:-1]
    mat._shape = (mat._shape[0]-1, mat._shape[1])


def trim_matrix(mat, i):
    """
    Trims a matrix by deleting both a row and a column in a matrix (warning: inefficient)

    :param mat: matrix
    :param i: row index
    :return:
    """
    mat = mat.copy()
    mat = mat.tocsr()
    delete_row_csr(mat, i)
    mat = mat.transpose()
    mat = mat.tocsr()
    delete_row_csr(mat, i)
    mat = mat.transpose()
    return mat


def sparse_abs(sparse_matrix: spmat.csc_matrix) -> spmat.csc_matrix:
    """
    Recovers an absolute value of a sparse matrix

    :param sparse_matrix: sparse matrix for which we want to recover the absolute.
    :return: absolute of that matrix
    """
    # with warnings.catch_warnings():
    #     warnings.filterwarnings("ignore", "Changing the sparsity structure")

    # print('sparse_abs method unroll')
    # print(type(sparse_matrix), sparse_matrix.dtype, sparse_matrix.count_nonzero())

    sign = sparse_matrix.sign()  # csc
    ret_mat = sparse_matrix.multiply(sign)  # csc

    return ret_mat


def build_sink_source_current_array(io_index_pair: Tuple[int, int],
                                    shape: Tuple[int, int]) -> spmat.csc_matrix:
    """
    converts index pair to a solver-compatible array

    :param shape: shape of the conductance matrix
    :param io_index_pair: pair of indexes where sinks/bioflow pair is
    :param splu: boolean flag e are using a version of splu
    """

    if switch_to_splu:
        io_array = np.zeros((shape[0]-1, 1))
        if io_index_pair[0] < io_index_pair[1]:
            io_array[io_index_pair[0], 0] = 1.0
        else:
            io_array[io_index_pair[0]-1, 0] = 1.0
        return spmat.csc_matrix(io_array)
    else:
        io_array = np.zeros((shape[0], 1))
        io_array[io_index_pair[0], 0], io_array[io_index_pair[1], 0] = (1.0, -1.0)
        return spmat.csc_matrix(io_array)


def get_potentials(conductivity_laplacian: spmat.csc_matrix,
                   io_index_pair: Tuple[int, int],
                   shared_solver: Union[chmd.Factor, None]) -> Union[chmd.Factor, np.array]:
    """
    Recovers voltages based on the conductivity Laplacian and the IO array

    :param conductivity_laplacian:
    :param io_index_pair:

    :return: array of potential in each node or solver - depending on method used
    """
    # technically, Cholesky is not the best solver. In approximation it is good enough and is faster
    i, j = io_index_pair

    if switch_to_splu:
        io_array = build_sink_source_current_array(io_index_pair, conductivity_laplacian.shape)
        local_conductivity_laplacian = trim_matrix(conductivity_laplacian, io_index_pair[1])
        local_conductivity_laplacian += spmat.eye(*local_conductivity_laplacian.shape) * line_loss
        # log.info('starting splu computation')
        solver = splu(local_conductivity_laplacian)  # solver sharing for splu is impossible
        # log.info('splu computation done')
        voltages = solver.solve(io_array.toarray())
        voltages = np.insert(voltages, j, 0, axis=0)

    else:
        io_array = build_sink_source_current_array(io_index_pair, conductivity_laplacian.shape)
        if not share_solver or shared_solver is None:
            # KNOWNBUG: unsolved issue with a thread/pool vs import vs cholesky problem here
            solver = chmd.cholesky(conductivity_laplacian, line_loss)
        else:
            solver = shared_solver
        voltages = solver(io_array)

    return voltages


def get_current_matrix(conductivity_laplacian: spmat.csc_matrix,
                       node_potentials: spmat.csc_matrix) -> Tuple[spmat.csc_matrix, spmat.csc_matrix]:
    """
    Recovers the current matrix based on the conductivity laplacian and voltages in each node.

    :param conductivity_laplacian:
    :param node_potentials:
    :return: matrix where M[i,j] = current intensity from i to j. Assymteric and Triangular
     superior iof the assymetric one. if current is from j to i, term is positive, otherwise
     it is negative.
    :rtype: scipy.sparse.csc_matrix
    """
    if switch_to_splu:
        diag_voltages = spmat.diags(node_potentials.T.tolist()[0], 0, format="csc")
    else:
        diag_voltages = spmat.diags(node_potentials.toarray().T.tolist()[0], 0, format="csc")

    # csc
    corr_conductance_matrix = conductivity_laplacian - \
                              spmat.diags(conductivity_laplacian.diagonal(), 0,
                                                           format="csc")

    _lft = diag_voltages.dot(corr_conductance_matrix)  # csc
    _rgt = corr_conductance_matrix.dot(diag_voltages)  # csc

    currents = _lft - _rgt  # csc

    abs_current = sparse_abs(currents)  # csc
    abs_current_t = abs_current.transpose()  # csr
    abs_current_t = abs_current_t.tocsc()  # csc
    alt_currents = abs_current + abs_current_t  # csc
    # delayed triangularization for performance reasons
    # tri_currents = spmat.triu(alt_currents, format='csc')  # csc
    tri_currents = alt_currents

    return alt_currents, tri_currents


def get_current_through_nodes(triu_current_matrix: spmat.csc_matrix) -> List[np.float64]:
    """
    Recovers current flowing through each node

    :param triu_current_matrix: non-redundant (i.e. triangular superior/lower) matrix of
    currents through a conduction system
    :return : current through the individual nodes based on the current matrix as defined in
     the get_current_matrix module
    :rtype: numpy.array
    :raise Exception: if the incoming vs outcoming current flow is different
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Changing the sparsity structure")

    triu_current_matrix = triu_current_matrix.tocsc()
    positive_current = spmat.csc_matrix(triu_current_matrix.shape)
    positive_current[triu_current_matrix > 0.0] = triu_current_matrix[triu_current_matrix > 0.0]
    # that's the part raising a warning about sparsity. not critical for performance given it's
    # outside the main loop

    incoming_current = np.array((positive_current + positive_current.T).sum(axis=1)).flatten() / 2

    if node_current_in_debug:
        negative_current = spmat.csc_matrix(triu_current_matrix.shape)
        negative_current[triu_current_matrix < 0.0] = triu_current_matrix[triu_current_matrix < 0.0]
        outgoing_current = np.array((negative_current + negative_current.T).sum(axis=1)).flatten() / 2

        if np.any(outgoing_current):
            log.critical('incoming', incoming_current)
            log.critical('outgoing', outgoing_current)
            log.critical('diff', incoming_current + outgoing_current)
            log.critical(np.any(incoming_current + outgoing_current > 0.))
            raise Exception('debug: assumption failed....')

    ret = list(incoming_current)

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
    diagonal_pad = spmat.diags(pad_array, 0, format="csc")
    re_laplacian = copy(laplacian)
    re_laplacian = diagonal_pad.dot(re_laplacian.dot(diagonal_pad))
    re_laplacian = re_laplacian - \
        spmat.diags(re_laplacian.diagonal(), 0, format="csc")
    d = (-re_laplacian.sum(axis=0)).tolist()[0]
    re_laplacian = re_laplacian + spmat.diags(d, 0, format="csc")

    return re_laplacian


def edge_current_iteration(conductivity_laplacian: spmat.csc_matrix,
                           index_pair: Tuple[int, int],
                           shared_solver: Union[chmd.Factor, None] = None,
                           reach_limiter=None) -> (np.float64, spmat.csc_matrix):
    """
    Master edge current retriever

    :param conductivity_laplacian:
    :param shared_solver:
    :param index_pair:
    :param reach_limiter:
    :return: potential_difference, triu_current
    """
    if shared_solver is None and switch_to_splu is False and reach_limiter is None:
        log.warning('edge current computation could be accelerated by using a shared solver')

    if reach_limiter:
        conductivity_laplacian = laplacian_reachable_filter(conductivity_laplacian, reach_limiter)
        shared_solver = None  # solver cannot be shared because the reach filter changes solver

    i, j = index_pair

    voltages = get_potentials(conductivity_laplacian, index_pair, shared_solver)

    # if solver:
    #     if switch_to_splu:
    #         pass  # because the solver was not passed to start with
    #     else:
    #         io_array = build_sink_source_current_array((i, j), conductivity_laplacian.shape)  # csc
    #         voltages = solver(io_array)  # csc;
    # else:
    #     voltages = get_potentials(conductivity_laplacian, (i, j))
    #     # problem: are we retrieving anything with i, j when we cut the laplacian?
    #     if switch_to_splu:
    #         voltages = np.insert(voltages, j, 0, axis=0)

    _, current = get_current_matrix(conductivity_laplacian, voltages)  # csc, csc;
    potential_diff = abs(voltages[i, 0] - voltages[j, 0])  # np.float64

    return potential_diff, current


def main_flow_calc_loop(conductivity_laplacian: np.array,
                        sample: List[Tuple[int, float]],
                        secondary_sample: Union[List[Tuple[int, float]], None] = None,
                        cancellation: bool = True, potential_dominated: bool = True,
                        sparse_rounds: int = -1,
                        memory_source=None,
                        potential_diffs_remembered: bool = False,
                        thread_hex: str = '______',
                        flow_calculation_method=general_flow):
    """
    master method for all the required edge current calculations

    :param conductivity_laplacian: conductivity laplacian
    :param sample: (index, weight) between which to calculate the flow
    :param secondary_sample: (index, weight) for star-like or biparty calculation
    :param cancellation: if the total current is normalized to the sampe number
    :param potential_dominated: if the total current is normalized to potential
    :param sparse_rounds: to which depth is the sampling performed. if < 1, none will be
    :param memory_source: if we are using memoization, source to look it
    :param potential_diffs_remembered: if the difference of potentials between nodes is remembered
    :param thread_hex: debugging id of the thread in which the sampling is going on
    :param flow_calculation_method: the function that converts the sample signature (sample,
        secondary_sample, sparse_rounds) into a list of ((index, weight), (index weight)) tuples
    :return:
    """

    if not memory_source_allowed:
        memory_source = None

    log.info('thread hex: %s; master edge current starting to sample with %s nodes; cancellation: '
             '%s;'
             ' potential-dominated %s; sparse_rounds %s'
             % (thread_hex, len(sample), cancellation,
                potential_dominated, sparse_rounds))

    # log.debug('debug: parameters supplied to main_flow_calc_loop: '
    #          'conductivity_laplacian: %s,\n'
    #          'sample: %s,\n'
    #          'cancellation: %s,\n'
    #          'potential_dominated: %s,\n'
    #          'sparse_sampling %s,\n'
    #          'sparse_sampling_depth %s,\n'
    #          'memory_source: %s,\n'
    #          'potential_diffs_remembered: %s,\n'
    #          'thread_hex: %s' % (
    #     type(conductivity_laplacian),
    #     type(sample),
    #     cancellation,
    #     potential_dominated,
    #     sparse_sampling,
    #     sparse_sampling_depth,
    #     type(memory_source),
    #     potential_diffs_remembered,
    #     thread_hex))

    # convert the arguments to proper structure:
    conductivity_laplacian = conductivity_laplacian.tocsc()

    # line_of_interest = 7493
    # log.info('mat_lines of interest max: %f, %f; %d, %d' %
    #          (np.max(conductivity_laplacian[line_of_interest, :]),
    #           np.max(conductivity_laplacian[:, line_of_interest]),
    #           np.sum((conductivity_laplacian[line_of_interest, :] != 0).astype(np.int)),
    #           np.sum((conductivity_laplacian[:, line_of_interest] != 0).astype(np.int))))

    # generate index list in agreement with the sparse_sampling strategy

    # if sparse_sampling:
    #     list_of_pairs = []
    #     for _ in repeat(None, sparse_sampling_depth):
    #         idx_list_c = copy(sample)
    #         random.shuffle(idx_list_c)
    #         list_of_pairs += list(zip(idx_list_c[:len(idx_list_c) // 2],
    #                              idx_list_c[len(idx_list_c) // 2:]))
    # else:
    #     list_of_pairs = [(i, j) for i, j in combinations(set(sample), 2)]

    list_of_pairs = flow_calculation_method(sample, secondary_sample, sparse_rounds)

    total_pairs = len(list_of_pairs)

    up_pair_2_voltage = {}
    current_accumulator = spmat.csc_matrix(conductivity_laplacian.shape)

    if share_solver and not switch_to_splu:
        importlib.reload(chmd)
        log.debug('Chmd reloaded')  # Correction tentative did not work.
        shared_solver = chmd.cholesky(conductivity_laplacian, line_loss)
    else:
        shared_solver = None


    # if not switch_to_splu and share_solver:
    #     cholesky_shared_solver = cholesky(conductivity_laplacian, line_loss)
    # else:
    #     cholesky_shared_solver = None

    # run the main loop on the list of indexes in agreement with the memoization strategy:
    breakpoints = 300
    previous_time = time()

    for counter, (i, j) in enumerate(list_of_pairs):

        mean_weight = (i[1] + j[1]) / 2.  # KNOWNBUG: not sure if it works if the weight of one is 0
        i, j = (i[0], j[0])

        if memory_source and tuple(sorted((i, j))) in list(memory_source.keys()):
            potential_diff, current_upper = memory_source[tuple(sorted((i, j)))]

        else:
            # types : np.float64, csc_matrix
            potential_diff, current_upper = \
                edge_current_iteration(conductivity_laplacian, (i, j), shared_solver=shared_solver)

        if potential_diffs_remembered:
            up_pair_2_voltage[tuple(sorted((i, j)))] = potential_diff

        # normalize to potential, if needed
        if potential_dominated:
            if potential_diff != 0:
                # csc_matrix
                current_upper = current_upper / potential_diff

            # warn if potential difference is null or close to it
            else:
                log.warning('pairwise flow. On indexes %s %s potential difference is null. %s',
                            i, j, 'Tension-normalization was aborted')

        # csc = csc + csc
        current_accumulator = current_accumulator + sparse_abs(current_upper) * mean_weight

        if counter % breakpoints == 0 and counter > 1:
            # TODO: [load bar]: the internal loop load bar goes here
            compops = float(breakpoints) / (time() - previous_time)
            mins_before_termination = (total_pairs - counter) / compops // 60
            finish_time = datetime.datetime.now() + datetime.timedelta(minutes=mins_before_termination)
            log.info("thread hex: %s; progress: %s/%s, current speed: %.2f compop/s, "
                     "time remaining: "
                     "%.0f "
                     "min, finishing: %s "
                     % (thread_hex, counter, total_pairs, compops, mins_before_termination,
                        finish_time.strftime("%m/%d/%Y, %H:%M:%S")))
            previous_time = time()

    current_accumulator = spmat.triu(current_accumulator)

    if cancellation:
        current_accumulator /= float(total_pairs)

    return current_accumulator, up_pair_2_voltage


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

    return inverter[1] / inverter[0], inverter[0]
