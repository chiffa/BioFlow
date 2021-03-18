"""
Linear Algebra routines used throughout the project
"""
from itertools import combinations_with_replacement, combinations, product
from random import randrange, shuffle

import matplotlib.pyplot as plt
import numpy as np

import scipy.sparse as spmat
from scipy.sparse import lil_matrix, triu
from scipy.sparse.linalg import eigsh

from sklearn.cluster import spectral_clustering
import warnings

from bioflow.utils.general_utils.useful_wrappers import time_it_wrapper
from bioflow.utils.log_behavior import get_logger
from bioflow.utils.dataviz import render_2d_matrix

log = get_logger(__name__)

# warnings.filterwarnings("ignore", message="Changing the sparsity structure")


# @time_it_wrapper
def normalize_laplacian(non_normalized_laplacian: spmat.lil_matrix) -> spmat.lil_matrix:
    """
    Performs a normalization of a laplacian of a graph

    :param non_normalized_laplacian: non-normalized laplacian
    :type non_normalized_laplacian: scipy.sparse matrix
    :return:
    """
    non_normalized_laplacian = non_normalized_laplacian.tolil()
    diagonal_terms = np.array(non_normalized_laplacian.diagonal().tolist())
    diagonal_terms[diagonal_terms == 0] = 1
    diagonal_terms = np.reshape(
        diagonal_terms, (non_normalized_laplacian.shape[0], 1))
    normalisation_matrix = np.sqrt(np.dot(diagonal_terms, diagonal_terms.T))
    matrix_to_return = spmat.lil_matrix(
        non_normalized_laplacian /
        normalisation_matrix)
    return matrix_to_return


def average_off_diag_in_sub_matrix(matrix, indexes, show=False):
    """
    Calculates the average of off-diagonal elements in a submatrix defined by the indexes

    :param matrix: matrix
    :param indexes: index set defining the submatrix
    :param show: if the results are rendered for a show
    :return:
    """
    extracted_matrix = lil_matrix((len(indexes), len(indexes)))
    main2local_idx = dict((j, i) for i, j in enumerate(sorted(indexes)))
    # inverts the sorting of the indexes?

    for index1, index2 in combinations_with_replacement(indexes, 2):
        extracted_matrix[
            main2local_idx[index1],
            main2local_idx[index2]] = matrix[
            index1,
            index2]
        extracted_matrix[
            main2local_idx[index2],
            main2local_idx[index1]] = matrix[
            index2,
            index1]

    if show:
        plt.imshow(
            extracted_matrix.toarray(),
            cmap='jet',
            interpolation="nearest")
        plt.colorbar()
        plt.show()

    if extracted_matrix.shape[0] > 1:
        retvalue = np.sum(np.triu(extracted_matrix.toarray(), 1)) / \
            (extracted_matrix.shape[0] * (extracted_matrix.shape[0] - 1) / 2)
        return retvalue
    else:
        return 0


def average_interset_linkage(matrix,
                             index_sets) -> np.float64:
    """
    Calculates the average linkage between several index sets within a matrix

    :param matrix: matrix of linkage
    :param index_sets: sets between which we want to calculate the linkage
    :return: average linkage between sets
    """
    accumulator = 0
    counter = 0
    for set1, set2 in combinations(index_sets, 2):
        for idx1, idx2 in product(set1, set2):
            accumulator += matrix[idx2, idx1]
            counter += 1

    return accumulator / counter


def show_eigenvals_and_eigenvects(eigenvals, eigenvects, biggest_limit,
                                  annotation, eigenvect_indexing_names=None):
    """
    Shows eigenvalues and eigenvectors supplied as parameters. If no names associated to specific
    eigenvectors positions are provided, will just print eigenvalues and exist.
    Else, will print statistics for the

    :param eigenvals: eigenvalues
    :param eigenvects: eigenvectors
    :param biggest_limit: how many biggest eigenvectors we would like to see
    :param annotation: additional annotation we would like to record in logs
    :param eigenvect_indexing_names: dict mapping indexes to names
    :return: None, this is a print-representation method only
    """
    log.info(annotation)
    log.info('########################')
    # if no index characters are supplied, just prints all the the eigenvalues
    # and exits
    if not eigenvect_indexing_names:
        for eigval in eigenvals.tolist():
            log.info(eigval)
        log.info('<============================>\n')
        return None

    absolute = np.absolute(eigenvects)
    biggest_eigenvects = np.argsort(absolute, axis=0)[-biggest_limit:, :]
    super_list = []
    # iterate through the  the eigenvalues, sums the contribution of each position in the N biggest
    # eigenvectors and adds it ot the log console
    for i, eigval in enumerate(reversed(eigenvals.tolist())):
        string_to_render = [
            'associated eigval: %s' %
            str(eigval),
            'Index \t value \t Descriptor']
        local_set = set()
        # the couple of lines below a kind a little bit esoteric...
        for index, value in reversed(list(zip(biggest_eigenvects[:, i],
                                         absolute[biggest_eigenvects[:, i], i]))):
            if int(value ** 2 * 1000) > 0:
                local_set.add(index)
                string_to_render.append('\t'.join([str(index), str(
                    '{0:.4f}'.format(value ** 2)), str(eigenvect_indexing_names[index])]))
        string_to_render.append('<==================>')
        # log the string we want to render only if it is the first time local set has been
        # encountered
        if local_set not in super_list:
            super_list.append(local_set)
            log.info('\n'.join(string_to_render))

    return None


# @time_it_wrapper
def analyze_eigenvects(
        non_normalized_laplacian,
        num_first_eigenvals_to_analyse,
        index_chars,
        permutations_limiter=1e7):
    """
    Performs an analysis of eigenvalues and eigenvectors of a non-normalized laplacian matrix
    in order to detect if anything found is statistically significant compared to a permutation
    control.

    :param non_normalized_laplacian:
    :param num_first_eigenvals_to_analyse:
    :param index_chars: dict mapping indexes to names
    :param permutations_limiter: maximal number of permutations
    :return:
    """
    log.debug('analyzing the laplacian with %s items and %s non-zero elements',
              non_normalized_laplacian.shape[0] ** 2,
              len(non_normalized_laplacian.nonzero()[0]))

    # normalize the laplacian
    normalized_laplacian = normalize_laplacian(non_normalized_laplacian)
    true_eigenvals, true_eigenvects = eigsh(
        normalized_laplacian, num_first_eigenvals_to_analyse)

    # extract non-zero indexes of matrix off-diagonal terms
    triangular_upper = lil_matrix(triu(normalized_laplacian))
    triangular_upper.setdiag(0)
    triangular_upper_nonzero_idxs = triangular_upper.nonzero()

    log.info("reassigning the indexes for %s items, with %s non-zero elements",
             triangular_upper.shape[0] ** 2,
             len(triangular_upper_nonzero_idxs[0]))

    # permute the off-diagonal terms
    nonzero_coordinates = list(zip(triangular_upper_nonzero_idxs[0].tolist(),
                              triangular_upper_nonzero_idxs[1].tolist()))
    shuffle(nonzero_coordinates)

    if len(nonzero_coordinates) > permutations_limiter:
        nonzero_coordinates = nonzero_coordinates[:permutations_limiter]

    for i, j in nonzero_coordinates:
        k = randrange(1, triangular_upper.shape[0] - 1)
        l = randrange(k + 1, triangular_upper.shape[0])
        triangular_upper[
            i, j], triangular_upper[
            k, l] = (
            triangular_upper[
                k, l], triangular_upper[
                    i, j])
    # rebuild normalized laplacian matrix from diagonal matrix
    off_diag_laplacian = triangular_upper + triangular_upper.T
    diag_laplacian = [- item for sublist in off_diag_laplacian.sum(
        axis=0).tolist() for item in sublist]
    off_diag_laplacian.setdiag(diag_laplacian)
    normalized_permutation_control = normalize_laplacian(off_diag_laplacian)
    # recompute the eigenvalues
    control_eigenvals, control_eigenvects = eigsh(
        normalized_permutation_control, num_first_eigenvals_to_analyse)
    show_eigenvals_and_eigenvects(true_eigenvals, true_eigenvects, 20,
                                  'true laplacian', index_chars)
    show_eigenvals_and_eigenvects(
        control_eigenvals,
        control_eigenvects,
        20,
        'random')


def cluster_nodes(dist_laplacian: spmat.lil_matrix,
                  clusters: int = 3,
                  show: bool = False) -> np.array:
    """
    Clusters indexes together based on a distance lapalacian. Basicall a wrapper for a spectral
    clustering from scipy package

    :param dist_laplacian: laplacian of internode distance
    :param clusters: desired number of clusters
    :param show: if we want to show the laplacian we are performing the computation on
    :return:
    """
    norm_laplacian = normalize_laplacian(dist_laplacian)
    norm_laplacian.setdiag(0)
    norm_laplacian = -norm_laplacian
    if show:
        render_2d_matrix(norm_laplacian.toarray(), "Node clustering normalized laplacian")

    labels = spectral_clustering(
        norm_laplacian,
        n_clusters=clusters,
        eigen_solver='arpack')

    return np.reshape(labels, (dist_laplacian.shape[0], 1))


if __name__ == "__main__":
    pass
