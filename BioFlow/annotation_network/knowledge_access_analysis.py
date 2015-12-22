"""
Set of methods responsible for knowledge analysis
"""
import pickle
from collections import namedtuple
from copy import copy
from csv import reader
from multiprocessing import Pool
from random import shuffle

import numpy as np
from matplotlib import pyplot as plt

from BioFlow.algorithms_bank.conduction_routines import perform_clustering
from BioFlow.annotation_network.BioKnowledgeInterface import GeneOntologyInterface
from BioFlow.main_configs import UP_rand_samp, Dumps, analysis_set_bulbs_ids,\
    background_set_bulbs_ids
from BioFlow.molecular_network.InteractomeInterface import InteractomeInterface
from BioFlow.utils.io_Routines import undump_object
from BioFlow.utils.linalg_routines import analyze_eigenvects
from BioFlow.utils.dataviz import kde_compute
from BioFlow.utils.log_behavior import get_logger

log = get_logger(__name__)


# TODO: refactor, this is currently a wrapper method used in the
def get_go_interface_instance():
    """
    Generates a Matrix_Knowledge_DB interface for the use in the spawner. If

    :return: a GO_interface object
    """
    ################################
    # Attention, manual switch here:
    ################################

    go_interface_instance = GeneOntologyInterface(_filter,
                                                  get_background(),
                                                  _correlation_factors,
                                                  True, 3)
    # go_interface_instance = GO_Interface(_filter, interactome_interface_instance.all_uniprots_bulbs_id_list, _correlation_factors, True, 3)
    go_interface_instance.load()
    log.info(go_interface_instance.pretty_time())
    return go_interface_instance


def spawn_sampler(sample_size_list_plus_iteration_list_plus_args):
    """
    Spawns a sampler initialized from the default GO_Interface

    :param sample_size_list_plus_iteration_list_plus_args: combined list of sample sizes
     and iterations (required for Pool.map usage)
    """
    go_interface_instance = sample_size_list_plus_iteration_list_plus_args[4]
    if go_interface_instance is None:
        go_interface_instance = get_go_interface_instance()

    sample_size_list = sample_size_list_plus_iteration_list_plus_args[0]
    iteration_list = sample_size_list_plus_iteration_list_plus_args[1]
    sparse_rounds = sample_size_list_plus_iteration_list_plus_args[2]
    chromosome_specific = sample_size_list_plus_iteration_list_plus_args[3]
    go_interface_instance.randomly_sample(
        sample_size_list,
        iteration_list,
        sparse_rounds,
        chromosome_specific)


def spawn_sampler_pool(
        pool_size,
        sample_size_list,
        iterations_list_per_pool,
        sparse_rounds=False,
        chromosome_specific=False,
        go_interface_instance=None):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param iterations_list_per_pool: number of iterations performing the pooling of the samples
     in each list
    :param sparse_rounds:
    :param chromosome_specific:
    :param go_interface_instance:
    """
    p = Pool(pool_size)
    payload = [
        (sample_size_list,
         iterations_list_per_pool,
         sparse_rounds,
         chromosome_specific,
         go_interface_instance)]
    p.map(spawn_sampler, payload * pool_size)


def select(tri_array, array_column, selection_span):
    """
    Convenient small function to select a from tri_array all the elements where the\ column
    number array_column is within the selection span

    :param tri_array: the matrix on which we will be performing the selection
    :param array_column: column on which the selection span will be applied
    :param selection_span: span for which we are going to keep the column.
    :return:
    """
    selector = np.logical_and(
        selection_span[0] < tri_array[
            array_column, :], tri_array[
            array_column, :] < selection_span[1])
    if not any(selector):
        return np.array([[0.0, 0.0, 0.0]])
    decvec = tri_array[:, selector]
    return decvec


def show_correlations(
        tri_corr_array,
        mean_correlations,
        eigenvalues,
        selector,
        test_tri_corr_array,
        test_mean_correlation,
        eigenvalue,
        re_samples,
        go_interface_instance=None):

    # TODO: there is a lot of repetition depending on which values are the biggest,
    # test-setted or real setted. In all, we should be able to reduce it to two functions:
    # scatterplot and histogram with two sets that should go into the dataviz module
    """
    A general function that performs demonstration of an example of random samples of the
     same size as our sample
    and of our sample and conducts the statistical tests on whether any of nodes or
     functional groups in our sample are non-random

    :param tri_corr_array: [[current, informativity, confusion_potential], ...] -
    characteristics of the random samples
    :param mean_correlations: [[cluster size, average internode connection], ...] -
    characteristics of clustering random samples with the same parameters
    :param eigenvalues: eigenvalues associated to the interconnection matrix of random samples
    :param selector: range on which we would like to visually zoom and plot a histogram
    :param test_tri_corr_array: [[current, informativity, confusion_potential], ...] -
    characteristics of the true sample. If none, nothing happens
    :param test_mean_correlation: [[cluster size, average internode connection], ...] -
    characteristics of clustering the true sample
    :param eigenvalue: eigenvalues associated to the interconnection matrix of the true sample
    :param re_samples: how many random samples we analyzed for the default model
    :param go_interface_instance:
    :return:
    """
    if go_interface_instance is None:
        go_interface_instance = get_go_interface_instance()

    inf_sel = (go_interface_instance.calculate_informativity(selector[0]),
               go_interface_instance.calculate_informativity(selector[1]))

    plt.figure()

    plt.subplot(331)
    plt.title('current through nodes')
    bins = np.linspace(tri_corr_array[0, :].min(),
                       tri_corr_array[0, :].max(),
                       100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(tri_corr_array[0, :].min(), test_tri_corr_array[0, :].min(
        )), max(tri_corr_array[0, :].max(), test_tri_corr_array[0, :].max()), 100)
    plt.hist(tri_corr_array[0, :],
             bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(test_tri_corr_array[0, :],
                 bins=bins, histtype='step', log=True, color='r')

    plt.subplot(332)
    plt.title('test current vs pure informativity')
    plt.scatter(tri_corr_array[1, :], tri_corr_array[0, :])
    if test_tri_corr_array is not None:
        plt.scatter(
            test_tri_corr_array[1, :],
            test_tri_corr_array[0, :],
            color='r', alpha=0.5)
    plt.axvspan(inf_sel[0], inf_sel[1], facecolor='0.5', alpha=0.3)

    plt.subplot(333)
    plt.title('test current v.s. confusion potential')
    plt.scatter(tri_corr_array[2, :], tri_corr_array[0, :])
    if test_tri_corr_array is not None:
        plt.scatter(
            test_tri_corr_array[2, :],
            test_tri_corr_array[0, :],
            color='r', alpha=0.5)
    plt.axvspan(selector[0], selector[1], facecolor='0.5', alpha=0.3)

    plt.subplot(334)
    plt.title('Gaussian KDE current_info')
    estimator_function = kde_compute(tri_corr_array[(1, 0), :], 50, re_samples)
    current_info_rel = None
    if test_tri_corr_array is not None:
        current_info_rel = estimator_function(test_tri_corr_array[(1, 0), :])

    plt.subplot(335)
    plt.title('GO_term pure informativity')
    bins = np.linspace(
        tri_corr_array[1, :].min(),
        tri_corr_array[1, :].max(),
        100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(tri_corr_array[1, :].min(),
                               test_tri_corr_array[1, :].min()),
                           max(tri_corr_array[1, :].max(),
                               test_tri_corr_array[1, :].max()),
                           100)
    plt.hist(tri_corr_array[1, :],
             bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(test_tri_corr_array[1, :],
                 bins=bins, histtype='step', log=True, color='r')

    plt.subplot(336)
    plt.title('Density of current in the highlighted area')
    bins = np.linspace(select(tri_corr_array, 2, selector)[0, :].min(),
                       select(tri_corr_array, 2, selector)[0, :].max(),
                       100)
    if test_tri_corr_array is not None:
        bins = np.linspace(
            min(select(tri_corr_array, 2, selector)[0, :].min(),
                select(test_tri_corr_array, 2, selector)[0, :].min()),
            max(select(tri_corr_array, 2, selector)[0, :].max(),
                select(test_tri_corr_array, 2, selector)[0, :].max()),
            100)

    plt.hist(select(tri_corr_array, 2, selector)[0, :],
             bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(select(test_tri_corr_array, 2, selector)[0, :],
                 bins=bins, histtype='step', log=True, color='r')

    # this property is better off viewed as a scatterplot of true points and
    # default points
    plt.subplot(337)
    plt.title('Clustering correlation')
    # plt.scatter(mean_correlations[0, :], mean_correlations[1, :], color = 'b')
    estimator_function = kde_compute(mean_correlations[(0, 1), :], 50, re_samples)
    cluster_props = None
    if test_mean_correlation is not None:
        plt.scatter(test_mean_correlation[0, :],
                    test_mean_correlation[1, :],
                    color='k', alpha=0.8)
        cluster_props = estimator_function(test_mean_correlation[(0, 1), :])

    plt.subplot(338)
    plt.title('Eigvals_hist')
    bins = np.linspace(eigenvalues.min(), eigenvalues.max(), 100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(eigenvalues.min(), eigenvalue.min()),
                           max(eigenvalues.max(), eigenvalue.max()),
                           100)
    plt.hist(eigenvalues, bins=bins, histtype='step', color='b')
    if eigenvalue is not None:
        plt.hist(eigenvalue.tolist() * 3, bins=bins, histtype='step', color='r')

    plt.subplot(339)
    plt.title('confusion potential')
    bins = np.linspace(tri_corr_array[2, :].min(),
                       tri_corr_array[2, :].max(),
                       100)
    if test_tri_corr_array is not None:
        bins = np.linspace(min(tri_corr_array[2, :].min(),
                               test_tri_corr_array[2, :].min()),
                           max(tri_corr_array[2, :].max(),
                               test_tri_corr_array[2, :].max()),
                           100)
    plt.hist(tri_corr_array[2, :],
             bins=bins, histtype='step', log=True, color='b')
    if test_tri_corr_array is not None:
        plt.hist(test_tri_corr_array[2, :],
                 bins=bins, histtype='step', log=True, color='r')

    plt.show()

    # pull the groups corresponding to non-random associations.
    return current_info_rel, cluster_props


def compare_to_blank(
        blank_model_size,
        zoom_range_selector,
        real_knowledge_interface=None,
        p_val=0.05,
        sparse_rounds=False,
        clusters=3,
        go_interface_instance=None):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param blank_model_size: the number of uniprots in the blanc model
    :param zoom_range_selector: tuple representing the coverage range for which we would want
     to see the histogram of current distributions
    :param real_knowledge_interface: The GO_Interface that has run the current computation
    :param p_val: desired p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used with
     the number of rounds equal to the number
    :param clusters: specifies the number of clusters we want to have
    :param go_interface_instance:
    :return: None if no significant nodes, the node and group characterisitc dictionaries
     otherwise
    """
    if go_interface_instance is None:
        go_interface_instance = get_go_interface_instance()

    md5_hash = go_interface_instance.md5_hash()

    curr_inf_conf_general = []
    count = 0
    mean_correlation_accumulator = []
    eigenvalues_accumulator = []

    log.info("requested: size: %s, sys_hash: %s, sparse_rounds: %s,",
             blank_model_size, md5_hash,  sparse_rounds)

    log.info("samples found to test against: \t %s",
             UP_rand_samp.find({'size': blank_model_size,
                                'sys_hash': md5_hash,
                                'sparse_rounds': sparse_rounds}).count())

    background_sample = UP_rand_samp.find(
            {'size': blank_model_size, 'sys_hash': md5_hash, 'sparse_rounds': sparse_rounds})

    for i, sample in enumerate(background_sample):
        _, node_currents = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])
        _, _, mean_correlations, eigenvalues = perform_clustering(tensions, 3, show=False)
        mean_correlation_accumulator.append(np.array(mean_correlations))
        eigenvalues_accumulator.append(eigenvalues)
        dict_system = go_interface_instance.format_node_props(node_currents)
        curr_inf_conf = list(dict_system.itervalues())
        curr_inf_conf_general.append(np.array(curr_inf_conf).T)
        count = i
        if dict_system == {}:
            del mean_correlation_accumulator[-1]
            del curr_inf_conf_general[-1]
            del eigenvalues_accumulator[-1]
            print "exceptional state: nothing in dict_system."

    # This part declares the pre-operators required for the verification of a
    # real sample

    final = np.concatenate(tuple(curr_inf_conf_general), axis=1)
    final_mean_correlations = np.concatenate(tuple(mean_correlation_accumulator), axis=0).T
    final_eigenvalues = np.concatenate(tuple(eigenvalues_accumulator), axis=0).T
    curr_inf_conf = None
    mean_correlations = None
    eigenvalue = None
    group2avg_off_diag = None
    go_node_ids = None
    dict_system = None
    if real_knowledge_interface:
        node_currents = real_knowledge_interface.node_current
        dict_system = go_interface_instance.format_node_props(node_currents)
        curr_inf_conf_tot = np.array(
            [[int(key)] + list(val) for key, val in dict_system.iteritems()]).T
        go_node_ids, curr_inf_conf = (
            curr_inf_conf_tot[
                0, :], curr_inf_conf_tot[
                (1, 2, 3), :])
        log.info('blank comparison: %s', curr_inf_conf.shape)
        group2avg_off_diag, _, mean_correlations, eigenvalue = perform_clustering(
            real_knowledge_interface.UP2UP_voltages, clusters)

    log.info('stats on %s samples', count)

    # TODO: We could and should separate the visualisation from the gaussian
    # estimators computation
    r_nodes, r_groups = show_correlations(
        final, final_mean_correlations, final_eigenvalues,
        zoom_range_selector, curr_inf_conf, mean_correlations.T, eigenvalue.T, count)

    go_node_char = namedtuple(
        'Node_Char', [
            'current', 'informativity', 'confusion_potential', 'p_value'])
    group_char = namedtuple(
        'Group_Char', [
            'UPs', 'num_UPs', 'average_connection', 'p_value'])
    if r_nodes is not None:
        not_random_nodes = [str(int(GO_id))
                            for GO_id in go_node_ids[r_nodes < p_val].tolist()]
        not_random_groups = np.concatenate(
            (group2avg_off_diag,
             np.reshape(r_groups, (3, 1))),
            axis=1)[r_groups < p_val].tolist()
        not_random_groups = [group_char(*nr_group)
                             for nr_group in not_random_groups]
        # basically the second element below are the nodes that contribute to the information
        #  flow through the node that is considered as non-random
        dct = dict((GO_id,
                    tuple([go_node_char(*(dict_system[GO_id] + r_nodes[go_node_ids ==
                                                                       float(GO_id)].tolist())),
                           list(set(go_interface_instance.GO2UP_Reachable_nodes[GO_id]
                                    ).intersection(set(real_knowledge_interface.analytic_Uniprots
                                                       )))]))
                   for GO_id in not_random_nodes)

        return sorted(dct.iteritems(), key=lambda x: x[
                      1][0][3]), not_random_groups

    return None, None, None


def decide_regeneration():
    """
    A script to decide at what point it is better to recompute a new a network rather
    then go through the time it requires to be upickled.
    The current decision is that for the samples of the size of ~ 100
    reached_uniprots_bulbs_id_list, we are better off unpickling from 4
    and more by factor 2 and by factor 10 from 9
    Previous experiments have shown that memoization with pickling incurred no noticeable
    delay on samples of up to
    50 UPs, but that the storage limit on mongo DB was rapidly exceeded, leading us to
    create an allocated dump file.
    """
    sample_root = []
    rooot_copy = copy(sample_root)
    go_interface_instance = get_go_interface_instance()
    go_interface_instance.set_uniprot_source(sample_root)
    go_interface_instance.build_extended_conduction_system()
    go_interface_instance.export_conduction_system()
    log.info('decide_regeneration 1: %s', go_interface_instance.pretty_time())
    for i in range(2, 9):
        shuffle(rooot_copy)
        go_interface_instance.export_subsystem(sample_root, rooot_copy[:i ** 2])
        log.info('decide_regeneration 2: %s, retrieve \t %s',
                 i ** 2, go_interface_instance.pretty_time())
        go_interface_instance.set_uniprot_source(rooot_copy[:i ** 2])
        go_interface_instance.build_extended_conduction_system(memoized=False)
        go_interface_instance.export_conduction_system()
        log.info('decide_regeneration 3: %s, redo: \t %s',
                 i ** 2, go_interface_instance.pretty_time())


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
        counter += sample_size * sample ** 2 / operations_per_sec
        log.info("Computing a sample of %s proteins would take %s secs",
                 sample, "{0:.2f}".format(sample ** 2 / operations_per_sec))
        log.info("Repeated %s times: %s h. Total time after phase: %s h",
                 sample_size,
                 "{0:.2f}".format(sample_size * sample ** 2 / operations_per_sec / 3600),
                 "{0:.2f}".format(counter / 3600))
    return counter


def linearly_independent_go_groups(size):
    """
    Performs the analysis of linearly independent GO groups

    :param size: number of most prominent eigenvectors contributing to an eigenvalue
    we would want to see
    """
    go_interface_instace = get_go_interface_instance()
    go_interface_instace.undump_independent_linear_sets()
    char_indexes = dict(
        (key,
         (len(
             go_interface_instace.GO2UP_Reachable_nodes[value]),
             go_interface_instace.GO_Legacy_IDs[value],
             go_interface_instace.GO_Names[value])) for key,
        value in go_interface_instace.Num2GO.iteritems())
    print go_interface_instace.pretty_time()
    analyze_eigenvects(go_interface_instace.Indep_Lapl, size, char_indexes)


def auto_analyze(source=None, go_interface_instance=None, processors=3, desired_depth=24):
    """
    Automatically analyzes the GO annotation of the RNA_seq results.

    :param source:
    :param go_interface_instance:
    :param processors:
    :param desired_depth:
    """
    if source is None:
        dumplist = undump_object(Dumps.RNA_seq_counts_compare)
    else:
        dumplist = [source]

    if desired_depth % processors != 0:
        desired_depth = desired_depth / processors + 1
    else:
        desired_depth = desired_depth / processors

    # noinspection PyTypeChecker
    for my_list in dumplist:
        if go_interface_instance is None:
            go_interface_instance = get_go_interface_instance()

        go_interface_instance.set_uniprot_source(my_list)
        log.info("auto_analyze_1: %s", len(go_interface_instance.analytic_uniprots))

        if len(go_interface_instance.analytic_uniprots) < 200:
            spawn_sampler_pool(processors,
                               [len(go_interface_instance.analytic_uniprots)],
                               [desired_depth],
                               go_interface_instance=go_interface_instance)
            go_interface_instance.build_extended_conduction_system()
            nr_nodes, nr_groups = compare_to_blank(
                len(go_interface_instance.analytic_uniprots),
                [1100, 1300],
                go_interface_instance,
                p_val=0.9, go_interface_instance=go_interface_instance)
        else:
            sampling_depth = max(200 ** 2 / len(go_interface_instance.analytic_uniprots), 5)

            log.info('lenght: %s \t sampling depth: %s \t, estimated_time: %s',
                     len(go_interface_instance.analytic_uniprots),
                     sampling_depth,
                     len(go_interface_instance.analytic_uniprots) * sampling_depth / 2 / 6 / 60)

            spawn_sampler_pool(processors,
                               [len(go_interface_instance.analytic_uniprots)],
                               [desired_depth],
                               sparse_rounds=sampling_depth,
                               go_interface_instance=go_interface_instance)

            go_interface_instance.build_extended_conduction_system(sparse_samples=sampling_depth)
            go_interface_instance.export_conduction_system()
            nr_nodes, nr_groups = compare_to_blank(
                len(go_interface_instance.analytic_uniprots),
                [1100, 1300],
                go_interface_instance,
                p_val=0.9, sparse_rounds=sampling_depth,
                go_interface_instance=go_interface_instance)

        go_interface_instance.export_conduction_system()

        for group in nr_groups:
            print group
        for node in nr_nodes:
            print node


def build_blank(length, depth, sparse_rounds=False):
    """
    Builds a blank set for the current analysis system

    :param length:
    :param depth:
    :param sparse_rounds:
    :return:
    """
    go_interface_instace = get_go_interface_instance()
    md5_hash = go_interface_instace.md5_hash()
    if UP_rand_samp.find({'size': length,
                          'sys_hash': md5_hash,
                          'sparse_rounds': sparse_rounds}).count() < depth:
        spawn_sampler_pool(4, [length], [depth])


def run_analysis(group):
    """
    Performs the whole analysis round retrieving

    :param group:
    :return:
    """
    go_interface_instance = get_go_interface_instance()
    go_interface_instance.set_uniprot_source(group)
    go_interface_instance.build_extended_conduction_system()
    go_interface_instance.export_conduction_system()
    nr_nodes, nr_groups = compare_to_blank(
        len(group), [1000, 1200], go_interface_instance, p_val=0.9)
    for group in nr_groups:
        print 'run analysis 1: %s' % group
    for node in nr_nodes:
        print 'run analysis 2: %s' % node


def get_source():
    """ retrieves bulbs ids for the elements for the analyzed group """
    source_bulbs_ids = []
    with open(analysis_set_bulbs_ids) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            source_bulbs_ids = source_bulbs_ids + row
    source_bulbs_ids = [ret for ret in source_bulbs_ids]
    return source_bulbs_ids


def get_background():
    """ retrieves bulbs ids for the elements in the background group """
    background_bulbs_ids = []
    with open(background_set_bulbs_ids) as src:
        csv_reader = reader(src)
        for row in csv_reader:
            background_bulbs_ids = background_bulbs_ids + row
    background_bulbs_ids = [ret for ret in background_bulbs_ids]
    return background_bulbs_ids


if __name__ == "__main__":

    _filter = ['biological_process']
    _correlation_factors = (1, 1)
    interactome_interface_instance = InteractomeInterface(True, False)
    interactome_interface_instance.fast_load()

    KG = GeneOntologyInterface(_filter,
                               interactome_interface_instance.all_uniprots_bulbs_id_list,
                               (1, 1), True, 3)
    KG.load()
    print KG.pretty_time()
    auto_analyze(get_source(), KG)
