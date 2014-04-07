__author__ = 'ank'

from scipy.sparse import lil_matrix, triu
from scipy.sparse.linalg import eigsh
import numpy as np
from copy import copy
from itertools import repeat
from random import randrange

def Lapl_normalize(non_normalized_Laplacian):
    """
    performs a normalization of a laplacian of a graph

    :param non_normalized_Laplacian: non-normalized laplacian
    :type non_normalized_Laplacian: scipy.sparse matrix
    :return:
    """
    normalized_Laplacian = lil_matrix(copy(non_normalized_Laplacian))
    diagterms = normalized_Laplacian.diagonal().tolist()
    for i, term in enumerate(diagterms):
        normalized_Laplacian[i, :] = normalized_Laplacian[i, :]/np.sqrt(term)
        normalized_Laplacian[:, i] = normalized_Laplacian[:, i]/np.sqrt(term)
    return normalized_Laplacian

def show_eigenvals_and_eigenvects(eigvals, eigvects, biggest_limit, annotation):
    print annotation, eigvals
    # absolute_vals = np.absolute(eigvals)
    # biggest_vals = np.argsort(absolute_vals, axis=0)[-biggest_limit:, :]
    # print biggest_vals
    # for eigval in biggest_vals
    #
    # norm = np.linalg.norm(MG.cond_eigenvects, axis = 1)
    # sbigest_args = np.argsort(norm)[- N_of_eigenvectors_to_characterise:]
    #
    # for arg in reversed(sbigest_args):
    #     print 'Matrix row: %s, \t Row Norm: %s, \t Description: %s' % (arg,  norm[arg],  MG.get_descriptor_for_index(arg))
    #     Retname.add(str(MG.get_descriptor_for_index(arg)[1]))


def analyze_eigvects(non_normalized_Laplacian, num_first_eigvals_to_analyse):

    # normalize the laplacian
    normalized_Laplacian = Lapl_normalize(non_normalized_Laplacian)
    # compute the eigenvalues and storre them
    true_eigenvals, true_eigenvects = eigsh(normalized_Laplacian, num_first_eigvals_to_analyse)
    # permute randomly the off-diagonal terms
    triag_u = lil_matrix(triu(non_normalized_Laplacian))
    triag_u.setdiag(0)
    for _ in repeat(None, non_normalized_Laplacian.shape[0]**2/2):
        i = randrange(1, non_normalized_Laplacian.shape[0]-1)
        j = randrange(i+1, non_normalized_Laplacian.shape[0])
        triag_u[i,j], triag_u[j,i] = (triag_u[j,i], triag_u[i,j])

    # recompute the diagonal terms
    fullmat = triag_u + triag_u.T
    diagterms = [item for sublist in fullmat.sum(axis=0).tolist() for item in sublist]
    fullmat.setdiag(diagterms)
    # recompute the normalized matrix
    normalized_rand = Lapl_normalize(fullmat)
    # recompute the eigenvalues
    rand_eigenvals, rand_eigenvects = eigsh(normalized_rand, num_first_eigvals_to_analyse)

    show_eigenvals_and_eigenvects(true_eigenvals, true_eigenvects, 10, 'true')
    show_eigenvals_and_eigenvects(rand_eigenvals, rand_eigenvects, 10, 'random')

    # normation and return of associated terms




if __name__ == "__main__":
    test_lapl = lil_matrix(np.zeros((3, 3)))
    test_lapl.setdiag([1, 2, 3])
    test_lapl[1, 2] = 2
    test_lapl[2, 1] = 2
    test_lapl[0, 2] = 1
    test_lapl[2, 0] = 1
    print test_lapl.toarray()
    print Lapl_normalize(test_lapl).toarray()


    analyze_eigvects(test_lapl,2)
