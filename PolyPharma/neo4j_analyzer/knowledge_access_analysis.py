"""
Builds on the knowledge
"""
__author__ = 'ank'

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
from PolyPharma.configs import UP_rand_samp
from PolyPharma.Utils.dataviz import kde_compute
from PolyPharma.Utils.Linalg_routines import analyze_eigvects
from PolyPharma.neo4j_analyzer.Conduction_routines import perform_clustering


filtr = ['biological_process']
corrfactors = (1, 1)
pprinter = PrettyPrinter(indent=4)
MG = MatrixGetter(True, False)
MG.fast_load()


def KG_gen():
    KG = GO_Interface(filtr, MG.Uniprot_complete, corrfactors, True, 3)
    KG.load()
    print KG.pretty_time()
    return KG


def spawn_sampler(sample_size_list_plus_iteration_list_plus_args):
    """
    Spawns a sampler initalized from the default GO_Interface

    :param sample_size_list_plus_iteration_list: combined list of sample swizes and iterations (requried for Pool.map usage)
    """
    KG = KG_gen()
    sample_size_list = sample_size_list_plus_iteration_list_plus_args[0]
    iteration_list = sample_size_list_plus_iteration_list_plus_args[1]
    sparse_rounds = sample_size_list_plus_iteration_list_plus_args[2]
    chromosome_specific = sample_size_list_plus_iteration_list_plus_args[3]
    KG.randomly_sample(sample_size_list, iteration_list, sparse_rounds, chromosome_specific)


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


def select(tri_array, array_column, selection_span):
    """
    Convinient small function to select a from tri_array all the elements where the column number array_column
    is within the selection span

    :param tri_array: the matrix on which we will be performing the selection
    :param array_column: column on which the selection span will be applied
    :param selection_span: span for which we are going to keep the column.
    :return:
    """
    selector = np.logical_and(selection_span[0]< tri_array[array_column, :], tri_array[array_column, :]< selection_span[1])
    if not any(selector):
        return np.array([[0.0,0.0,0.0]])
    decvec = tri_array[:, selector]
    return decvec


def show_corrs(tri_corr_array, meancorrs, eigvals, selector, test_tri_corr_array, test_meancorr, eigval, resamples):

    # TODO: there is a lot of repetition depending on which values are the biggest, test-setted or real setted.
    # In all, we should be able to reduce it to two functions: scatterplot and histogram with two sets that
    # should go into the dataviz module

    """
    A general function that performs demonstration of an example of random samples of the same size as our sample
    and of our sample and conducts the statistical tests on wherther any of nodes or functional groups in our
    sample are non-random

    :param tri_corr_array: [[current, informativity, confusion_potential], ...] - characteristics of the random samples
    :param meancorrs: [[cluster size, average internode connection], ...] - characteristics of clustering random samples with the same parameters
    :param eigvals: eigenvalues associated to the interconnection matrix of random samples
    :param selector: range on which we would like to visually zoom and plot a histogram
    :param test_tri_corr_array: [[current, informativity, confusion_potential], ...] - characteristics of the true sample. If none, nothing happens
    :param test_meancorr: [[cluster size, average internode connection], ...] - characteristics of clustering the true sample
    :param eigval: eigenvalues associated to the interconnection matrix of the true sample
    :param resamples: how many random samples we analyzed for the default model
    :return:
    """
    KG = KG_gen()
    inf_sel = (KG._infcalc(selector[0]), KG._infcalc(selector[1]))

    plt.figure()

    plt.subplot(331)
    plt.title('current through nodes')
    bins = np.linspace(tri_corr_array[0, :].min(), tri_corr_array[0, :].max(), 100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(tri_corr_array[0, :].min(), test_tri_corr_array[0, :].min()),
                           max(tri_corr_array[0, :].max(), test_tri_corr_array[0, :].max()),
                           100)
    plt.hist(tri_corr_array[0, :], bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(test_tri_corr_array[0, :], bins=bins, histtype='step', log=True, color='r')

    plt.subplot(332)
    plt.title('test current vs pure informativity')
    plt.scatter(tri_corr_array[1, :], tri_corr_array[0, :])
    if test_tri_corr_array is not None:
        plt.scatter(test_tri_corr_array[1, :], test_tri_corr_array[0, :], color='r' , alpha=0.5)
    plt.axvspan( inf_sel[0], inf_sel[1], facecolor='0.5', alpha=0.3)

    plt.subplot(333)
    plt.title('test current v.s. confusion potential')
    plt.scatter(tri_corr_array[2, :], tri_corr_array[0, :])
    if test_tri_corr_array is not None:
        plt.scatter(test_tri_corr_array[2, :], test_tri_corr_array[0, :], color='r', alpha=0.5)
    plt.axvspan( selector[0], selector[1], facecolor='0.5', alpha=0.3)

    plt.subplot(334)
    plt.title('Gaussian KDE current_info')
    estimator_function = kde_compute(tri_corr_array[(1,0), :], 50, resamples)
    current_info_rel = None
    if test_tri_corr_array is not None:
        current_info_rel = estimator_function(test_tri_corr_array[(1,0),:])

    plt.subplot(335)
    plt.title('GO_term pure informativity')
    bins = np.linspace(tri_corr_array[1, :].min(), tri_corr_array[1, :].max(), 100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(tri_corr_array[1, :].min(), test_tri_corr_array[1, :].min()),
                           max(tri_corr_array[1, :].max(), test_tri_corr_array[1, :].max()),
                           100)
    plt.hist(tri_corr_array[1, :], bins=100, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(test_tri_corr_array[1, :], bins=100, histtype='step', log=True, color='r')

    plt.subplot(336)
    plt.title('Density of current in the highlighted area')
    bins = np.linspace(select(tri_corr_array, 2, selector)[0, :].min(),
                       select(tri_corr_array, 2, selector)[0, :].max(),
                       100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(select(tri_corr_array, 2, selector)[0, :].min(),
                               select(test_tri_corr_array, 2, selector)[0, :].min()),
                           max(select(tri_corr_array, 2, selector)[0, :].max(),
                               select(test_tri_corr_array, 2, selector)[0, :].max()),
                           100)
    plt.hist(select(tri_corr_array, 2, selector)[0, :], bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(select(test_tri_corr_array, 2, selector)[0, :], bins=bins, histtype='step', log=True, color='r')

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
    if test_tri_corr_array is not None:
        bins = np.linspace(min(eigvals.min(), eigval.min()),
                           max(eigvals.max(), eigval.max()),
                           100)
    plt.hist(eigvals, bins=bins, histtype='step', color = 'b')
    if eigval is not None:
        plt.hist(eigval.tolist()*3, bins=bins, histtype='step', color = 'r')

    plt.subplot(339)
    plt.title('confusion potential')
    bins = np.linspace(tri_corr_array[2, :].min(), tri_corr_array[2, :].max(), 100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(tri_corr_array[2, :].min(), test_tri_corr_array[2, :].min()),
                           max(tri_corr_array[2, :].max(), test_tri_corr_array[2, :].max()),
                           100)
    plt.hist(tri_corr_array[2, :], bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(test_tri_corr_array[2, :], bins=bins, histtype='step', log=True, color='r')

    plt.show()

    # pull the groups corresponding to non-random associations.
    return current_info_rel, cluster_props


def compare_to_blanc(blanc_model_size, zoom_range_selector, real_knowledge_interface = None, p_val=0.05, sparse_rounds=False):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param blanc_model_size:
    :param zoom_range_selector:
    :param real_knowledge_interface:
    :param p_val:
    :param sparse_rounds:
    :type sparse_rounds: int
    :return:
    """
    KG = KG_gen()
    MD5_hash = KG._MD5hash()

    curr_inf_conf_general = []
    count = 0
    meancorr_acccumulator = []
    eigval_accumulator = []

    # this part computes the items required for the creation of a blanc model
    print "requested", 'size:', blanc_model_size, 'sys_hash:', MD5_hash, 'sparse_rounds:', sparse_rounds
    print "samples found to test against:\t", UP_rand_samp.find({'size': blanc_model_size, 'sys_hash' : MD5_hash, 'sparse_rounds':sparse_rounds}).count()
    for i, sample in enumerate(UP_rand_samp.find({'size': blanc_model_size, 'sys_hash' : MD5_hash, 'sparse_rounds':sparse_rounds})):
        _, node_currs = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])
        _, _, meancorr, eigvals = perform_clustering(tensions, 3, show=False)
        meancorr_acccumulator.append(np.array(meancorr))
        eigval_accumulator.append(eigvals)
        Dic_system = KG.format_Node_props(node_currs)
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
    GO_node_ids = None
    Dic_system = None
    if real_knowledge_interface:
        node_currs = real_knowledge_interface.node_current
        Dic_system = KG.format_Node_props(node_currs)
        curr_inf_conf_tot = np.array([[int(key)]+list(val) for key, val in Dic_system.iteritems()]).T
        GO_node_ids, curr_inf_conf = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1,2,3), :])
        group2avg_offdiag, _, meancorr, eigval = perform_clustering(real_knowledge_interface.UP2UP_voltages, 3)


    print "stats on %s samples" % count

    # TODO: We could and should separate the visualisation from the gaussian estimators computation
    r_nodes, r_groups = show_corrs(final, final_meancorrs, final_eigvals, zoom_range_selector, curr_inf_conf, meancorr.T, eigval.T, count)


    GO_node_char = namedtuple('Node_Char',['current','informativity', 'confusion_potential', 'p_value'])
    Group_char = namedtuple('Group_Char', ['UPs','num_UPs','average_connection','p_value'])
    if r_nodes is not None:
        not_random_nodes = [str(int(GO_id)) for GO_id in GO_node_ids[r_nodes < p_val].tolist()]
        not_random_groups = np.concatenate((group2avg_offdiag, np.reshape(r_groups, (3,1))), axis=1)[r_groups < p_val].tolist()
        not_random_groups = [ Group_char(*(nr_group)) for nr_group in not_random_groups]
        # basically the second element below are the nodes that contribute to the information flow through the node that is considered as
        # non-random
        dct = dict((GO_id,
                      tuple([GO_node_char(*(Dic_system[GO_id] + r_nodes[GO_node_ids == float(GO_id)].tolist())),
                              list(set(KG.GO2UP_Reachable_nodes[GO_id]).intersection(set(real_knowledge_interface.analytic_Uniprots)))]))
                     for GO_id in not_random_nodes)

        return  sorted(dct.iteritems(), key=lambda x:x[1][0][3]), not_random_groups

    return  None


def decide_regeneration():
    """
    A script to decide at what point it is better to recompute a new a network rather then go through the time it
    requires to be upickled.
    The current decision is that for the samples of the size of ~ 100 Uniprots, we are better off unpickling from 4
    and more by factor 2 and by factor 10 from 9
    Previous experimets have shown that memoization with pickling incurred no noticeable delay on samples of up to
    50 UPs, but that the storage limit on mongo DB was rapidly exceeded, leading us to create an allocated dump file.
    """

    sample_root = ['530239', '921394', '821224', '876133', '537471', '147771', '765141', '783757', '161100', '996641',
    '568357', '832606', '888857', '443125', '674855', '703465', '770454', '585061', '767349', '454684', '476323',
    '890779', '699374', '699085', '926841', '719433', '979188', '750252', '884148', '452226', '510869', '934804',
    '450711', '654463', '475017', '836128', '869961', '833908', '748293', '642129', '511971', '450103', '465344',
    '664249', '759667', '479667', '945097', '934005', '474459', '616764', '993605', '151251', '881579', '1010120',
    '567103', '177132', '914246', '818797', '734031', '983957', '988876', '907270', '764944', '457147', '574367',
    '605567', '950635', '184544', '372652', '440372', '630427', '446382', '195073', '790029', '941318', '572041',
    '469852', '1009569', '969215', '571784', '794977', '991385', '511515', '592947', '517667', '746802', '685187',
    '877373', '860329', '589231', '1013595', '679330', '808630', '774665', '663924', '615588', '497135', '628832',
    '841054', '657304']
    rooot_copy = copy(sample_root)
    KG = KG_gen()
    KG.set_Uniprot_source(sample_root)
    KG.build_extended_conduction_system()
    KG.export_conduction_system()
    print KG.pretty_time()
    for i in range(2, 9):
        shuffle(rooot_copy)
        KG.export_subsystem(sample_root,rooot_copy[:i**2])
        print i**2, 'retrieve: \t', KG.pretty_time()
        KG.set_Uniprot_source(rooot_copy[:i**2])
        KG.build_extended_conduction_system(memoized = False)
        KG.export_conduction_system()
        print i**2, '    redo: \t', KG.pretty_time()


def get_estimated_time(samples, sample_sizes, operations_per_sec=2.2):
    """
    Short time to estimate the time required for the generation of random saples

    :param samples: size of the sample
    :param sample_sizes: times each sample size is sampled
    :param operations_per_sec: complex operations per second
    :return:
    """
    counter = 0
    for sample, sample_size in zip(samples, sample_sizes):
        counter += sample_size*sample**2/operations_per_sec
        print "Computing a sample of %s proteins would take %s secs. \t Repeated %s times: %s h. Total time after phase: %s h" %(
            sample, "{0:.2f}".format(sample**2/operations_per_sec),
            sample_size, "{0:.2f}".format(sample_size*sample**2/operations_per_sec/3600),
            "{0:.2f}".format(counter/3600)
        )
    return counter


def linindep_GO_groups(size):
    """
    Performs the analysis of linerly independent GO groups

    :param size:
    """
    KG = KG_gen()
    KG.undump_Indep_Linset()
    char_indexes = dict( (key, (len(KG.GO2UP_Reachable_nodes[value]), KG.GO_Legacy_IDs[value], KG.GO_Names[value])) for key, value in KG.Num2GO.iteritems())
    print KG.pretty_time()
    analyze_eigvects(KG.Indep_Lapl, size, char_indexes)

# TODO: write chromosome comparator


if __name__ == "__main__":
    # spawn_sampler(([10, 100], [2, 1]))
    # spawn_sampler_pool(4, [201], [10], sparse_rounds=100)

    # get_estimated_time([10, 25, 50, 100,], [15, 10, 10, 8,])
    # test_list = ['147875', '130437', '186024', '100154', '140777', '100951', '107645', '154772']
    # print len(test_list)
    KG = KG_gen()


    KG.randomly_sample([201],[1], chromosome_specific=15, sparse_rounds=100, memoized=True, No_add=True)
    print KG.current_accumulator.shape
    KG.export_conduction_system()
    nr_nodes, nr_groups = compare_to_blanc(201, [1000, 1200], KG, p_val=0.9, sparse_rounds=100)
    for group in nr_groups:
        print group
    for node in nr_nodes:
        print node


    # KG.set_Uniprot_source(test_list)
    # KG.build_extended_conduction_system()
    # KG.export_conduction_system()
    # nr_nodes, nr_groups = compare_to_blanc(len(test_list), [1000, 1200], KG, p_val=0.9)
    # for group in nr_groups:
    #     print group
    # for node in nr_nodes:
    #     print node

    # linindep_GO_groups(50)

    pass