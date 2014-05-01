__author__ = 'ank'


# reimplements the routines from the knowledge access analysis
# for the reactome elements

import pickle
import numpy as np

from copy import copy
from random import shuffle
from collections import namedtuple
from multiprocessing import Pool
from pprint import PrettyPrinter
from matplotlib import pyplot as plt
from Matrix_Interactome_DB_interface import  MatrixGetter
from Matrix_Knowledge_DB_Interface import GO_Interface
from PolyPharma.configs import Interactome_rand_samp
from PolyPharma.Utils.dataviz import kde_compute
from PolyPharma.Utils.Linalg_routines import analyze_eigvects
from PolyPharma.neo4j_analyzer.Conduction_routines import perform_clustering


pprinter = PrettyPrinter(indent=4)
MG = MatrixGetter(True, False)
MG.fast_load()


def MG_gen():
    MG = MatrixGetter(True, False)
    MG.fast_load()
    print MG.pretty_time()
    return MG

def spawn_sampler(sample_size_list_plus_iteration_list_plus_args):
    """
    Spawns a sampler initalized from the default GO_Interface

    :param sample_size_list_plus_iteration_list: combined list of sample swizes and iterations (requried for Pool.map usage)
    """
    MG = MG_gen()
    sample_size_list = sample_size_list_plus_iteration_list_plus_args[0]
    iteration_list = sample_size_list_plus_iteration_list_plus_args[1]
    sparse_rounds = sample_size_list_plus_iteration_list_plus_args[2]
    chromosome_specific = sample_size_list_plus_iteration_list_plus_args[3]
    MG.randomly_sample(sample_size_list, iteration_list, sparse_rounds, chromosome_specific)


def spawn_sampler_pool(pool_size, sample_size_list, interation_list_per_pool, sparse_rounds=False, chromosome_specific=False):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param interation_list_per_pool: number of iterations performing the pooling of the samples in each list
    :param sparse_rounds:
    :type sparse_rounds: int
    :param chromosome_specific:
    :type chromosome_specific: int
    """
    p = Pool(pool_size)
    payload = [(sample_size_list, interation_list_per_pool, sparse_rounds, chromosome_specific)]
    p.map(spawn_sampler, payload * pool_size)


def select(bi_array, array_column, selection_span):
    """
    Convinient small function to select a from tri_array all the elements where the column number array_column
    is within the selection span

    :param bi_array: the matrix on which we will be performing the selection
    :param array_column: column on which the selection span will be applied
    :param selection_span: span for which we are going to keep the column.
    :return:
    """
    selector = np.logical_and(selection_span[0]< bi_array[array_column, :], bi_array[array_column, :]< selection_span[1])
    if not any(selector):
        return np.array([[0.0,0.0,0.0]])
    decvec = bi_array[:, selector]
    return decvec


def show_corrs(bi_corr_array, meancorrs, eigvals, selector, test_bi_corr_array, test_meancorr, eigval, resamples):

    # TODO: there is a lot of repetition depending on which values are the biggest, test-setted or real setted.
    # In all, we should be able to reduce it to two functions: scatterplot and histogram with two sets that
    # should go into the dataviz module

    """
    A general function that performs demonstration of an example of random samples of the same size as our sample
    and of our sample and conducts the statistical tests on wherther any of nodes or functional groups in our
    sample are non-random

    :param bi_corr_array: [[current, informativity, confusion_potential], ...] - characteristics of the random samples
    :param meancorrs: [[cluster size, average internode connection], ...] - characteristics of clustering random samples with the same parameters
    :param eigvals: eigenvalues associated to the interconnection matrix of random samples
    :param selector: range on which we would like to visually zoom and plot a histogram
    :param test_bi_corr_array: [[current, informativity, confusion_potential], ...] - characteristics of the true sample. If none, nothing happens
    :param test_meancorr: [[cluster size, average internode connection], ...] - characteristics of clustering the true sample
    :param eigval: eigenvalues associated to the interconnection matrix of the true sample
    :param resamples: how many random samples we analyzed for the default model
    :return:
    """
    plt.figure()

    plt.subplot(331)
    plt.title('current through nodes')
    bins = np.linspace(bi_corr_array[0, :].min(), bi_corr_array[0, :].max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(bi_corr_array[0, :].min(), test_bi_corr_array[0, :].min()),
                           max(bi_corr_array[0, :].max(), test_bi_corr_array[0, :].max()),
                           100)
    plt.hist(bi_corr_array[0, :], bins=bins, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(test_bi_corr_array[0, :], bins=bins, histtype='step', log=True, color='r')

    plt.subplot(332)
    plt.title('test current vs degree')
    plt.scatter(bi_corr_array[1, :], bi_corr_array[0, :])
    if test_bi_corr_array is not None:
        plt.scatter(test_bi_corr_array[1, :], test_bi_corr_array[0, :], color='r' , alpha=0.5)
    plt.axvspan( selector[0], selector[1], facecolor='0.5', alpha=0.3)

    plt.subplot(333)
    plt.title('Currently empty')

    plt.subplot(334)
    plt.title('Gaussian KDE current_info')
    estimator_function = kde_compute(bi_corr_array[(1,0), :], 50, resamples)
    current_info_rel = None
    if test_bi_corr_array is not None:
        current_info_rel = estimator_function(test_bi_corr_array[(1,0),:])

    plt.subplot(335)
    plt.title('Node degree distribution')
    bins = np.linspace(bi_corr_array[1, :].min(), bi_corr_array[1, :].max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(bi_corr_array[1, :].min(), test_bi_corr_array[1, :].min()),
                           max(bi_corr_array[1, :].max(), test_bi_corr_array[1, :].max()),
                           100)
    plt.hist(bi_corr_array[1, :], bins=100, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(test_bi_corr_array[1, :], bins=100, histtype='step', log=True, color='r')

    plt.subplot(336)
    plt.title('Density of current in the highlighted area')
    if test_bi_corr_array is not None:
        bins = np.linspace(min(select(bi_corr_array, 1, selector)[0, :].min(),
                               select(test_bi_corr_array, 1, selector)[0, :].min()),
                           max(select(bi_corr_array, 1, selector)[0, :].max(),
                               select(test_bi_corr_array, 1, selector)[0, :].max()),
                           100)
    plt.hist(select(bi_corr_array, 1, selector)[0, :], bins=bins, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(select(test_bi_corr_array, 1, selector)[0, :], bins=100, histtype='step', log=True, color='r')

    # this property is better off viewed as a scatterplot of true points and default points
    plt.subplot(337)
    plt.title('Clustering correlation')
    # plt.scatter(meancorrs[0, :], meancorrs[1, :], color = 'b')
    estimator_function = kde_compute(meancorrs[(0,1), :], 50, resamples)
    cluster_props =None
    if test_meancorr is not None:
        plt.scatter(test_meancorr[0,:], test_meancorr[1,:], color = 'k', alpha=0.8)
        cluster_props = estimator_function(test_meancorr[(0,1),:])

    plt.subplot(338)
    plt.title('Eigvals_hist')
    bins = np.linspace(eigvals.min(), eigvals.max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(eigvals.min(), eigval.min()),
                           max(eigvals.max(), eigval.max()),
                           100)
    plt.hist(eigvals, bins=bins, histtype='step', color = 'b')
    if eigval is not None:
        plt.hist(eigval.tolist()*3, bins=bins, histtype='step', color = 'r')

    plt.subplot(339)
    plt.title('Currently empty')


    plt.show()

    # pull the groups corresponding to non-random associations.
    return current_info_rel, cluster_props


def compare_to_blanc(blanc_model_size, zoom_range_selector, real_interactome_interface = None, p_val=0.05, sparse_rounds=False, clusters=3):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :return:
    """
    MD5_hash = MG._MD5hash()

    curr_inf_conf_general = []
    count = 0
    meancorr_acccumulator = []
    eigval_accumulator = []

    print "samples found to test against:\t", Interactome_rand_samp.find({'size': blanc_model_size, 'sys_hash' : MD5_hash, 'sparse_rounds':sparse_rounds}).count()
    # this part computes the items required for the creation of a blanc model
    for i, sample in enumerate(Interactome_rand_samp.find({'size': blanc_model_size, 'sys_hash' : MD5_hash, 'sparse_rounds':sparse_rounds})):
        if sparse_rounds:
            raise Exception('blanc done on sparse rounds, clustering likely to be hazardous. Process interrupted')
        _, node_currs = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])
        _, _, meancorr, eigvals = perform_clustering(tensions, clusters, show=False)
        meancorr_acccumulator.append(np.array(meancorr))
        eigval_accumulator.append(eigvals)
        Dic_system = MG.format_Node_props(node_currs)
        curr_inf_conf = list(Dic_system.itervalues())
        curr_inf_conf_general.append(np.array(curr_inf_conf).T)
        count = i

    # This part declares the pre-operators required for the verification of a real sample

    final = np.concatenate(tuple(curr_inf_conf_general), axis=1)
    final_meancorrs = np.concatenate(tuple(meancorr_acccumulator), axis=0).T
    final_eigvals = np.concatenate(tuple(eigval_accumulator), axis=0).T
    curr_inf_conf = None
    meancorr = None
    eigval = None
    group2avg_offdiag = None
    node_ids = None
    Dic_system = None
    if real_interactome_interface:
        node_currs = real_interactome_interface.node_current
        Dic_system = MG.format_Node_props(node_currs)
        curr_inf_conf_tot = np.array([[int(key)] + list(val) for key, val in Dic_system.iteritems()]).T
        node_ids, curr_inf_conf = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1,2), :])
        group2avg_offdiag, _, meancorr, eigval = perform_clustering(real_interactome_interface.UP2UP_voltages, clusters)

    print "stats on %s samples" % count


    # TODO: We could and should separate the visualisation from the gaussian estimators computation
    r_nodes, r_groups = show_corrs(final, final_meancorrs, final_eigvals, zoom_range_selector, curr_inf_conf, meancorr.T, eigval.T, count)


    Interactome_node_char = namedtuple('Node_Char',['name', 'current', 'degree', 'p_value'])
    Group_char = namedtuple('Group_Char', ['UPs','num_UPs','average_connection','p_value'])

    if r_nodes is not None:
        not_random_nodes = [str(int(node_id)) for node_id in node_ids[r_nodes < p_val].tolist()]
        not_random_groups = np.concatenate((group2avg_offdiag, np.reshape(r_groups, (3,1))), axis=1)[r_groups < p_val].tolist()
        not_random_groups = [ Group_char(*(nr_group)) for nr_group in not_random_groups]
        # basically the second element below are the nodes that contribute to the information flow through the node that is considered as
        # non-random
        dct = dict((nr_node_id, Interactome_node_char(MG.ID2displayName[nr_node_id],*(Dic_system[nr_node_id] + r_nodes[node_ids == float(nr_node_id)].tolist())))
                     for nr_node_id in not_random_nodes)

        return  sorted(dct.iteritems(), key=lambda x:x[1][3]), not_random_groups

    return  None


if __name__ == "__main__":
    spawn_sampler_pool(4, [150], [10])

    MG1 = MG_gen()

    MG1.randomly_sample([150], [1], chromosome_specific=15, No_add=True)
    nr_nodes, nr_groups = compare_to_blanc(150, [0.5, 0.6], MG1, p_val=0.9)
    MG1.export_conduction_system()
    for group in nr_groups:
        print group
    for node in nr_nodes:
        print node

    pass