"""
New analytical routines for the interactome
"""
import pickle
from collections import namedtuple
from csv import reader
from multiprocessing import Pool

import numpy as np
from matplotlib import pyplot as plt

from bioflow.algorithms_bank.conduction_routines import perform_clustering
from bioflow.main_configs import interactome_rand_samp, Outputs, Dumps, estimated_comp_ops
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.log_behavior import get_logger
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids

log = get_logger(__name__)


# TODO: factor that into the "retrieve" routine of the laplacian wrapper
def get_interactome_interface():
    """
    Retrieves an "InteractomeInterface" object

    :return:
    """
    interactome_interface_instance = InteractomeInterface(main_connex_only=True,
                                                          full_impact=False)
    interactome_interface_instance.fast_load()
    log.debug("get_interactome state e_p_u_b_i length: %s",
              len(interactome_interface_instance.entry_point_uniprots_bulbs_ids))
    log.info("interactome interface loaded in %s" % interactome_interface_instance.pretty_time())
    # is the case now
    return interactome_interface_instance


# TODO: defined only as an auxilary to the spawn_sampler_pool => inline that in the sampler pool
def spawn_sampler(sample_size_list_plus_iteration_list_plus_args):
    """
    Spawns a sampler initialized from the default GO_Interface

    :param sample_size_list_plus_iteration_list_plus_args: combined list of sample sizes and
    iterations (required for Pool.map usage)
    """
    interactome_interface_instance_arg = sample_size_list_plus_iteration_list_plus_args[4]
    if interactome_interface_instance_arg is None:
        interactome_interface_instance = get_interactome_interface()
    else:
        interactome_interface_instance = interactome_interface_instance_arg

    sample_size_list = sample_size_list_plus_iteration_list_plus_args[0]
    iteration_list = sample_size_list_plus_iteration_list_plus_args[1]
    sparse_rounds = sample_size_list_plus_iteration_list_plus_args[2]
    chromosome_specific = sample_size_list_plus_iteration_list_plus_args[3]  # TODO: remove that
    interactome_interface_instance.randomly_sample(
        sample_size_list,
        iteration_list,
        sparse_rounds,
        chromosome_specific)
    # TODO: remove chromosome-specificity. We are performing this analysis otherwise.


def spawn_sampler_pool(
        pool_size,
        sample_size_list,
        interaction_list_per_pool,
        sparse_rounds=False,
        chromosome_specific=False,
        interactome_interface_instance=None):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param interaction_list_per_pool: number of iterations performing the pooling of the samples
     in each list
    :param sparse_rounds:
    :param chromosome_specific:
    :param interactome_interface_instance:
    """
    # TODO: a better estimation of execution time
    process_pool = Pool(pool_size)
    payload = [
        (sample_size_list,
         interaction_list_per_pool,
         sparse_rounds,
         chromosome_specific,
         interactome_interface_instance)]
    log.debug('spawning the sampler with payload %s', payload)
    process_pool.map(spawn_sampler, payload * pool_size)


def local_indexed_select(bi_array, array_column, selection_span):
    """
    Convenient small function to select a from tri_array all the elements where the column
    number array_column is within the selection span

    :param bi_array: the matrix on which we will be performing the selection
    :param array_column: column number on which the selection span will be applied
    :param selection_span: span for which we are going to keep the column.
    """
    selector = np.logical_and(
        selection_span[0] < bi_array[array_column, :],
        bi_array[array_column, :] < selection_span[1])
    if not any(selector):
        return np.array([[0.0, 0.0, 0.0]])
    filtered_bi_array = bi_array[:, selector]
    return filtered_bi_array


# TODO: there is a lot of repetition depending on which values are the biggest,
# test-sets or real sets. In all, we should be able to reduce it to two functions:
# scatter plot and histogram with two sets that should go into the dataviz module
def show_test_statistics(
        bi_corr_array,
        mean_correlations,  # shutdown if sparse
        eigenvalues,        # shutdown if sparse
        selector,
        test_bi_corr_array,
        test_mean_corr,     # shutdown if sparse
        eigenvalue,         # shutdown fi sparse
        re_samples,
        sparse=False):
    """
    A general function that performs demonstration of an example of random samples of
     the same size as our sample and of our sample and conducts the statistical tests
     on wherther any of nodes or functional groups in our sample are non-random

    :param bi_corr_array: [[current, informativity, confusion_potential], ...] -
    characteristics of the random samples
    :param mean_correlations: [[cluster size, average internode connection], ...] -
    characteristics of clustering random samples with the same parameters
    :param eigenvalues: eigenvalues associated to the interconnection matrix of random samples
    :param selector: range on which we would like to visually zoom and plot a histogram
    :param test_bi_corr_array: [[current, informativity, confusion_potential], ...] -
    characteristics of the true sample. If none, nothing happens
    :param test_mean_corr: [[cluster size, average internode connection], ...] -
    characteristics of clustering the true sample
    :param eigenvalue: eigenvalues associated to the interconnection matrix of the true sample
    :param re_samples: how many random samples we analyzed for the default model
    :param sparse: True if we are showing test statistics of a sparse kernel run
    :return:
    """
    fig = plt.figure()
    fig.set_size_inches(30, 20)

    plt.subplot(331)
    plt.title('current through nodes')
    bins = np.linspace(
        bi_corr_array[0, :].min(),
        bi_corr_array[0, :].max(), 100)

    if test_bi_corr_array is not None:
        bins = np.linspace(min(bi_corr_array[0, :].min(), test_bi_corr_array[0, :].min()),
                           max(bi_corr_array[0, :].max(), test_bi_corr_array[0, :].max()),
                           100)
    plt.hist(bi_corr_array[0, :],
             bins=bins, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(test_bi_corr_array[0, :],
                 bins=bins, histtype='step', log=True, color='r')

    plt.subplot(332)
    plt.title('test current vs degree')
    plt.scatter(bi_corr_array[1, :], bi_corr_array[0, :])
    if test_bi_corr_array is not None:
        plt.scatter(test_bi_corr_array[1, :], test_bi_corr_array[0, :],
                    color='r', alpha=0.5)
    plt.axvspan(selector[0], selector[1], facecolor='0.5', alpha=0.3)

    plt.subplot(333)
    plt.title('Currently empty')

    plt.subplot(334)
    plt.title('Gaussian KDE current_info')
    estimator_function = kde_compute(bi_corr_array[(1, 0), :], 50, re_samples)
    current_info_rel = None
    if test_bi_corr_array is not None:
        current_info_rel = estimator_function(test_bi_corr_array[(1, 0), :])

    plt.subplot(335)
    plt.title('Node degree distribution')
    bins = np.linspace(bi_corr_array[1, :].min(), bi_corr_array[1, :].max(), 100)
    if test_bi_corr_array is not None:
        bins = np.linspace(min(bi_corr_array[1, :].min(), test_bi_corr_array[1, :].min()),
                           max(bi_corr_array[1, :].max(), test_bi_corr_array[1, :].max()),
                           100)
    plt.hist(bi_corr_array[1, :],
             bins=100, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(test_bi_corr_array[1, :],
                 bins=100, histtype='step', log=True, color='r')

    plt.subplot(336)
    plt.title('Density of current in the highlighted area')
    if test_bi_corr_array is not None:
        bins = np.linspace(
            min(local_indexed_select(bi_corr_array, 1, selector)[0, :].min(),
                local_indexed_select(test_bi_corr_array, 1, selector)[0, :].min()),
            max(local_indexed_select(bi_corr_array, 1, selector)[0, :].max(),
                local_indexed_select(test_bi_corr_array, 1, selector)[0, :].max()),
            100)
    plt.hist(local_indexed_select(bi_corr_array, 1, selector)[0, :],
             bins=bins, histtype='step', log=True, color='b')
    if test_bi_corr_array is not None:
        plt.hist(local_indexed_select(test_bi_corr_array, 1, selector)[0, :],
                 bins=100, histtype='step', log=True, color='r')

    # this property is better off viewed as a scatterplot of true points and
    # default points

    cluster_props = None

    plt.subplot(337)
    plt.title('Clustering correlation')
    if not sparse:
        # plt.scatter(mean_correlations[0, :], mean_correlations[1, :], color = 'b')
        estimator_function = kde_compute(mean_correlations[(0, 1), :], 50, re_samples)
        cluster_props = None
        if test_mean_corr is not None:
            plt.scatter(test_mean_corr[0, :], test_mean_corr[1, :],
                        color='k', alpha=0.8)
            cluster_props = estimator_function(test_mean_corr[(0, 1), :])

    plt.subplot(338)
    plt.title('Eigvals_hist')
    if not sparse:
        bins = np.linspace(eigenvalues.min(), eigenvalues.max(), 100)
        if test_bi_corr_array is not None:
            bins = np.linspace(min(eigenvalues.min(), eigenvalue.min()),
                               max(eigenvalues.max(), eigenvalue.max()),
                               100)
        plt.hist(eigenvalues, bins=bins, histtype='step', color='b')
        if eigenvalue is not None:
            plt.hist(eigenvalue.tolist() * 3, bins=bins, histtype='step', color='r')

    plt.subplot(339)
    plt.title('Currently empty')

    # plt.show()
    plt.savefig(Outputs.interactome_network_stats)

    # pull the groups corresponding to non-random associations.
    return current_info_rel, cluster_props


# TODO: correct the error: remove the clustering on the sparse round
def compare_to_blank(
        blank_model_size,
        zoom_range_selector,
        p_val=0.05,
        sparse_rounds=False,
        cluster_no=3,
        interactome_interface_instance=None):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param blank_model_size: the number of uniprots in the blank model
    :param zoom_range_selector: tuple representing the coverage range for which we would
     want to see the histogram of current distributions
    :param p_val: desired p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used
     with the number of rounds equal the integer value of that argument
    :param cluster_no: specifies the number of cluster_no we want to have
    :param interactome_interface_instance:
    :return: None if no significant nodes, the node and group characteristic
     dictionaries otherwise
    """
    if interactome_interface_instance is None:
        interactome_interface_instance = InteractomeInterface(True, False)
        interactome_interface_instance.fast_load()

    md5_hash = interactome_interface_instance.md5_hash()

    curr_inf_conf_general = []
    count = 0
    mean_correlation_accumulator = []
    eigenvalues_accumulator = []

    log.info("samples found to test against:\t %s" %
             interactome_rand_samp.find({'size': blank_model_size, 'sys_hash': md5_hash,
                                        'sparse_rounds': sparse_rounds}).count())

    # this part computes the items required for the creation of a blank model

    for i, sample in enumerate(interactome_rand_samp.find(
            {'size': blank_model_size, 'sys_hash': md5_hash, 'sparse_rounds': sparse_rounds})):
        if sparse_rounds:
            log.warning('Blank done on sparse rounds. Clustering will not be performed')
        _, node_currents = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])
        if not sparse_rounds:
            _, _, mean_correlations, eigvals = perform_clustering(
                tensions, cluster_no, show=False)
        else:
            mean_correlations = np.array([[(0, ), 0, 0]]*cluster_no)
            eigvals = np.array([-1]*cluster_no)
        mean_correlation_accumulator.append(np.array(mean_correlations))
        eigenvalues_accumulator.append(eigvals)
        dictionary_system = interactome_interface_instance.format_node_props(node_currents)
        curr_inf_conf = list(dictionary_system.itervalues())
        curr_inf_conf_general.append(np.array(curr_inf_conf).T)
        count = i

    # This part declares the pre-operators required for the verification of a
    # real sample

    final = np.concatenate(tuple(curr_inf_conf_general), axis=1)
    final_mean_correlations = np.concatenate(tuple(mean_correlation_accumulator), axis=0).T
    final_eigenvalues = np.concatenate(tuple(eigenvalues_accumulator), axis=0).T

    node_currents = interactome_interface_instance.node_current
    dictionary_system = interactome_interface_instance.format_node_props(node_currents)
    curr_inf_conf_tot = np.array(
        [[int(key)] + list(val) for key, val in dictionary_system.iteritems()]).T
    node_ids, curr_inf_conf = (curr_inf_conf_tot[0, :],
                               curr_inf_conf_tot[(1, 2), :])

    if not sparse_rounds:
        group2avg_offdiag, _, mean_correlations, eigenvalue = perform_clustering(
            interactome_interface_instance.UP2UP_voltages, cluster_no, 'Interactome clustering')

    else:
        group2avg_offdiag = np.array([[(0, ), 0, 0]]*cluster_no)
        mean_correlations = np.array([[0, 0]]*cluster_no)
        eigenvalue = np.array([-1]*cluster_no)

    log.info("stats on  %s samples" % count)

    # TODO: We could and should separate the visualisation from the gaussian
    # estimators computation
    r_nodes, r_groups = show_test_statistics(
        final, final_mean_correlations, final_eigenvalues,
        zoom_range_selector, curr_inf_conf, mean_correlations.T,
        eigenvalue.T, count, sparse_rounds)

    group_char = namedtuple(
        'Group_Char', [
            'UPs', 'num_UPs', 'average_connection', 'p_value'])

    if r_nodes is not None:
        not_random_nodes = [node_id for node_id in node_ids[r_nodes < p_val].tolist()]

        if not sparse_rounds:
            not_random_groups = np.concatenate(
                (group2avg_offdiag,
                 np.reshape(r_groups, (3, 1))), axis=1)[r_groups < p_val].tolist()
            not_random_groups = [group_char(*nr_group)
                                 for nr_group in not_random_groups]

        else:
            not_random_groups = []

        # basically the second element below are the nodes that contribute to the
        #  information flow through the node that is considered as non-random

        log.debug('debug, not random nodes: %s', not_random_nodes)
        log.debug('debug bulbs_id_disp_name: %s',
                  interactome_interface_instance.bulbs_id_2_display_name.items()[:10])

        node_char_list = [
            [int(nr_node_id),
             interactome_interface_instance.bulbs_id_2_display_name[nr_node_id]] +
            dictionary_system[nr_node_id] + r_nodes[node_ids == float(nr_node_id)].tolist()
            for nr_node_id in not_random_nodes]

        return sorted(node_char_list, key=lambda x: x[4]), not_random_groups

    return None, None


# TODO: check if background list works as expected
def auto_analyze(source_list,
                 desired_depth=24, processors=4,
                 background_list=None, skip_sampling=False):
    """
    Automatically analyzes the itneractome synergetic action of the RNA_seq results

    :param source_list::
    :param desired_depth:
    :param processors:
    :param background_list
    :param skip_sampling: if true, will skip background sampling step
    """
    # noinspection PyTypeChecker
    if desired_depth % processors != 0:
        desired_depth = desired_depth / processors + 1
    else:
        desired_depth = desired_depth / processors

    for _list in source_list:
        log.info('Auto analyzing list of interest: %s', len(_list))
        interactome_interface = get_interactome_interface()
        log.debug("retrieved interactome_interface instance e_p_u_b_i length: %s",
                  len(interactome_interface.entry_point_uniprots_bulbs_ids))
        interactome_interface.set_uniprot_source(list(_list))
        log.debug(" e_p_u_b_i length after UP_source was set: %s",
                  len(interactome_interface.entry_point_uniprots_bulbs_ids))

        interactome_interface.background = background_list
        if not skip_sampling:
            log.info("spawning a sampler for %s proteins @ %s compops/sec",
                     len(interactome_interface.entry_point_uniprots_bulbs_ids), estimated_comp_ops)

        if len(interactome_interface.entry_point_uniprots_bulbs_ids) < 60:

            if not skip_sampling:
                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(interactome_interface.entry_point_uniprots_bulbs_ids),
                         'full',
                         len(interactome_interface.entry_point_uniprots_bulbs_ids)**2 /
                         estimated_comp_ops / 60)

                spawn_sampler_pool(
                    processors,
                    [len(interactome_interface.entry_point_uniprots_bulbs_ids)],
                    [desired_depth],
                    interactome_interface_instance=None)

            interactome_interface.build_extended_conduction_system()
            nr_nodes, nr_groups = compare_to_blank(
                len(interactome_interface.entry_point_uniprots_bulbs_ids),
                [0.5, 0.6],
                p_val=0.9, interactome_interface_instance=interactome_interface)


        else:
            ceiling = min(205, len(interactome_interface.entry_point_uniprots_bulbs_ids))
            sampling_depth = max((ceiling-5) ** 2 /
                                 len(interactome_interface.entry_point_uniprots_bulbs_ids),
                                 5)

            if not skip_sampling:
                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(interactome_interface.entry_point_uniprots_bulbs_ids),
                         sampling_depth,
                         len(interactome_interface.entry_point_uniprots_bulbs_ids) *
                         sampling_depth / 2 / 60)

                spawn_sampler_pool(processors,
                                   [len(interactome_interface.entry_point_uniprots_bulbs_ids)],
                                   [desired_depth],
                                   sparse_rounds=sampling_depth,
                                   interactome_interface_instance=None)

            interactome_interface.build_extended_conduction_system(sparse_samples=sampling_depth)
            # interactome_interface.export_conduction_system()
            nr_nodes, nr_groups = compare_to_blank(
                len(interactome_interface.entry_point_uniprots_bulbs_ids),
                [0.5, 0.6],
                p_val=0.9, sparse_rounds=sampling_depth,
                interactome_interface_instance=interactome_interface)

        interactome_interface.export_conduction_system()
        for group in nr_groups:
            log.info(group)
        for node in nr_nodes:
            log.info('\t %s \t %s \t %s \t %s \t %s', *node)


if __name__ == "__main__":

    # pprinter = PrettyPrinter(indent=4)
    # interactome_interface_instance = MatrixGetter(True, False)
    # interactome_interface_instance.fast_load()

    # dumplist = undump_object(Dumps.RNA_seq_counts_compare)

    # MG1.randomly_sample([150], [1], chromosome_specific=15, No_add=True)
    # nr_nodes, nr_groups = compare_to_blanc(150, [0.5, 0.6], MG1, p_val=0.9)
    # MG1.export_conduction_system()
    # for group in nr_groups:
    #     print group
    # for node in nr_nodes:
    #     print node

    # source = get_source_bulbs_ids()
    # background_list = get_background_bulbs_ids()
    # auto_analyze([source], desired_depth=5, processors=6,
    #              background_list=background_list, skip_sampling=True)


    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=False)
    local_matrix.fast_load()
    #
    # spawn_sampler_pool(3, [2], [150], interactome_interface_instance=None)

    local_matrix.randomly_sample([195], [10], sparse_rounds=195)
