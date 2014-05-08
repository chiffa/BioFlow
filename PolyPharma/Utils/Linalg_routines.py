__author__ = 'ank'

from scipy.sparse import lil_matrix, triu
from scipy.sparse.linalg import eigsh
import numpy as np
from copy import copy
from itertools import combinations_with_replacement, combinations, product
from random import randrange, shuffle
from time import time
import matplotlib.pyplot as plt
from sklearn.cluster import spectral_clustering



def Lapl_normalize(non_normalized_Laplacian):
    """
    performs a normalization of a laplacian of a graph

    :param non_normalized_Laplacian: non-normalized laplacian
    :type non_normalized_Laplacian: scipy.sparse matrix
    :return:
    """
    normalized_Laplacian = lil_matrix(copy(non_normalized_Laplacian))
    diagterms = np.array(normalized_Laplacian.diagonal().tolist())
    diagterms[diagterms == 0] = 1
    diagterms = np.reshape(diagterms,(normalized_Laplacian.shape[0],1))
    normmat = np.sqrt(np.dot(diagterms, diagterms.T))
    retmat = lil_matrix(normalized_Laplacian/normmat)
    # retmat[retmat > 1.0] = 1.0
    # retmat[retmat <- 1.0] = -1.0
    return retmat


def submatrix(matrix, indexes, show = False):
    re_matrix = lil_matrix((len(indexes), len(indexes)))
    main2local_idx = dict( (j, i) for i, j in enumerate(sorted(indexes)))

    for index1, index2 in combinations_with_replacement(indexes, 2):
        re_matrix[main2local_idx[index1], main2local_idx[index2]] = matrix[index1, index2]
        re_matrix[main2local_idx[index2], main2local_idx[index1]] = matrix[index2, index1]

    if show:
        plt.imshow(re_matrix.toarray(), cmap='jet', interpolation="nearest")
        plt.colorbar()
        plt.show()

    if re_matrix.shape[0]>1:
        return np.sum(np.triu(re_matrix.toarray(), 1)) / (re_matrix.shape[0]*(re_matrix.shape[0]-1)/2)
    else:
        return 0

def remaineder_matrix(matrix, index_sets):
    accumulator = 0
    counter = 0
    for set1, set2 in combinations(index_sets,2):
        for idx1, idx2 in product(set1, set2):
            accumulator += matrix[idx2, idx1]
            counter+=1

    return accumulator/counter


def show_eigenvals_and_eigenvects(eigvals, eigvects, biggest_limit, annotation, idx_chars=None):
    print annotation
    print '########################'
    if not idx_chars:
        for eigval in eigvals.tolist():
            print eigval
        print '<============================>\n'
        return None

    absolute = np.absolute(eigvects)
    eigbiggest = np.argsort(absolute, axis=0)[-biggest_limit:, :]
    setmemory = []
    for i, eigval in enumerate(reversed(eigvals.tolist())):
        renderstring = []
        renderstring.append('associated eigval:'+str(eigval))
        renderstring.append('Index \t value \t Descriptor')
        mset = set()
        for index, value in reversed(zip(eigbiggest[:, i], absolute[eigbiggest[:, i], i])):
            if int(value**2*1000) > 0:
                mset.add(index)
                renderstring.append('\t'+str(index)+'\t'+str('{0:.4f}'.format(value**2))+'\t'+str(idx_chars[index]))
        renderstring.append('<==================>')
        if mset not in setmemory:
            setmemory.append(mset)
            print '\n'.join(renderstring)

    return None

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


def analyze_eigvects(non_normalized_Laplacian, num_first_eigvals_to_analyse, index_chars, permutations_limiter = 10000000, fudge=10e-10):
    # normalize the laplacian
    print 'analyzing the laplacian with %s items and %s non-zero elts' % (non_normalized_Laplacian.shape[0]**2, len(non_normalized_Laplacian.nonzero()[0]))
    t = time()
    init = time()
    normalized_Laplacian = Lapl_normalize(non_normalized_Laplacian)
    print time()-t
    t = time()
    # compute the eigenvalues and storre them
    true_eigenvals, true_eigenvects = eigsh(normalized_Laplacian, num_first_eigvals_to_analyse)
    print time()-t
    t = time()
    # permute randomly the off-diagonal terms
    triag_u = lil_matrix(triu(normalized_Laplacian))
    triag_u.setdiag(0)
    tnz = triag_u.nonzero()
    print "reassigning the indexes for %s items, with %s non-zero elts" % (triag_u.shape[0]**2, len(tnz[0]))
    eltsuite = zip(tnz[0].tolist(), tnz[1].tolist())
    shuffle(eltsuite)
    if eltsuite > permutations_limiter:
        eltsuite = eltsuite[:permutations_limiter]   # pb: we want it to affect any random number with reinsertion
    print time()-t
    t = time()
    # take a nonzero pair of indexes
    for i,j in eltsuite:
        # select randomly a pair of indexes and permute it
        k = randrange(1, triag_u.shape[0]-1)
        l = randrange(k+1, triag_u.shape[0])
        triag_u[i, j], triag_u[k,l] = (triag_u[k, l], triag_u[i, j])
    print time()-t
    t = time()
    # recompute the diagonal terms
    fullmat = triag_u + triag_u.T
    diagterms = [-item for sublist in fullmat.sum(axis=0).tolist() for item in sublist]
    fullmat.setdiag(diagterms)
    print time()-t
    t = time()
    # recompute the normalized matrix
    normalized_rand = Lapl_normalize(fullmat)
    # recompute the eigenvalues
    rand_eigenvals, rand_eigenvects = eigsh(normalized_rand, num_first_eigvals_to_analyse)
    print time()-t
    t = time()
    show_eigenvals_and_eigenvects(true_eigenvals, true_eigenvects, 20, 'true laplacian', index_chars)
    show_eigenvals_and_eigenvects(rand_eigenvals, rand_eigenvects, 20, 'random')
    print "final", time()-t, time()-init


def view_laplacian_off_terms(non_normalized_Laplacian):
    normalized_Laplacian = Lapl_normalize(non_normalized_Laplacian)
    triag_u = lil_matrix(triu(normalized_Laplacian))
    triag_u.setdiag(0)
    pre_arr = -triag_u[triag_u.nonzero()].toarray().flatten()
    arr = np.log10(pre_arr)
    plt.hist(arr, bins=100, log=True, histtype='step')
    plt.show()


def cluster_nodes(dist_laplacian, clusters=3, show=False):
    norm_laplacian = Lapl_normalize(dist_laplacian)
    norm_laplacian.setdiag(0)
    norm_laplacian = -norm_laplacian
    if show:
        plt.imshow(norm_laplacian.toarray(), cmap='jet', interpolation="nearest")
        plt.colorbar()
        plt.show()
    labels = spectral_clustering(norm_laplacian, n_clusters=clusters, eigen_solver='arpack')
    return np.reshape(labels, (dist_laplacian.shape[0], 1))



if __name__ == "__main__":
    test_lapl = lil_matrix(np.zeros((4, 4)))
    test_lapl.setdiag([1, 2, 3, 0])
    test_lapl[1, 2] = -2
    test_lapl[2, 1] = -2
    test_lapl[0, 2] = -1
    test_lapl[2, 0] = -1
    print test_lapl.toarray()
    print Lapl_normalize(test_lapl).toarray()

    analyze_eigvects(test_lapl, 2, {0:'zero', 1:'one', 2:'two', 3: 'three'})

    test_lapl[2, 3] = -0.5
    test_lapl[3, 2] = -0.5

    print cluster_nodes(test_lapl)
