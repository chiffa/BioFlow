"""
Set of methods responsible for knowledge analysis
"""
import pickle
from collections import namedtuple
from csv import reader
from multiprocessing import Pool
import traceback
import os

import numpy as np
from matplotlib import pyplot as plt
from csv import writer as csv_writer

from bioflow.algorithms_bank.conduction_routines import perform_clustering
from bioflow.annotation_network.BioKnowledgeInterface import GeneOntologyInterface
from bioflow.main_configs import annotome_rand_samp, Dumps, Outputs, estimated_comp_ops
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.io_routines import undump_object, get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.log_behavior import get_logger

log = get_logger(__name__)

_filter = ['biological_process']
_correlation_factors = (1, 1)
ref_param_set = tuple([_filter, [], (1, 1), True, 3])


# TODO: move ref_param_set to configs
def get_go_interface_instance(param_set=ref_param_set):
    """
    Generates a Matrix_Knowledge_DB interface for the use in the spawner. If

    :return: a GO_interface object
    """
    ################################
    # Attention, manual switch here:
    ################################

    go_interface_instance = GeneOntologyInterface(*param_set)
    # go_interface_instance = GO_Interface(_filter, interactome_interface_instance.
    # all_uniprots_neo4j_id_list, _correlation_factors, True, 3)
    go_interface_instance.load()
    log.info(go_interface_instance.pretty_time())
    return go_interface_instance


def spawn_sampler(sample_size_list_plus_iteration_list_plus_args):
    """
    Spawns a sampler initialized from the default GO_Interface

    :param sample_size_list_plus_iteration_list_plus_args: combined list of sample sizes
     and iterations (required for Pool.map usage)
    """
    try:
        go_interface_instance = sample_size_list_plus_iteration_list_plus_args[4]
        param_set = sample_size_list_plus_iteration_list_plus_args[5]
        if go_interface_instance is None:
            go_interface_instance = get_go_interface_instance(param_set)

        sample_size_list = sample_size_list_plus_iteration_list_plus_args[0]
        iteration_list = sample_size_list_plus_iteration_list_plus_args[1]
        sparse_rounds = sample_size_list_plus_iteration_list_plus_args[2]
        chromosome_specific = sample_size_list_plus_iteration_list_plus_args[3]
        go_interface_instance.randomly_sample(
            sample_size_list,
            iteration_list,
            sparse_rounds,
            chromosome_specific)

    except Exception as e:
        msg = "{}\n\nOriginal {}".format(e, traceback.format_exc())
        raise type(e)(msg)

    # TODO: choromosome specific has nothing to do here; should be managed otherwise


def spawn_sampler_pool(
        pool_size,
        sample_size_list,
        iterations_list_per_pool,
        sparse_rounds=False,
        chromosome_specific=False,
        go_interface_instance=None,
        param_set=ref_param_set):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param iterations_list_per_pool: number of iterations performing the pooling of the samples
     in each list
    :param sparse_rounds:
    :param chromosome_specific:
    :param go_interface_instance:
    :param param_set: set of parameters configuring the knowledge interface object
    """
    p = Pool(pool_size)
    payload = [
        (sample_size_list,
         iterations_list_per_pool,
         sparse_rounds,
         chromosome_specific,
         go_interface_instance,
         param_set)]
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
        go_interface_instance=None,
        sparse=False,
        param_set=ref_param_set):

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
    :param sparse:
    :param param_set:
    :return:
    """
    if go_interface_instance is None:
        go_interface_instance = get_go_interface_instance(param_set)

    inf_sel = (go_interface_instance.calculate_informativity(selector[0]),
               go_interface_instance.calculate_informativity(selector[1]))

    fig = plt.figure()
    fig.set_size_inches(30, 20)

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

    cluster_props = None

    plt.subplot(337)
    plt.title('Clustering correlation')
    if not sparse:
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
    if not sparse:
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

    # # plt.show()
    # plt.savefig(Outputs.knowledge_network_stats)

    # pull the groups corresponding to non-random associations.
    return current_info_rel, cluster_props


def compare_to_blank(
        blank_model_size,
        zoom_range_selector,
        p_val=0.05,
        sparse_rounds=False,
        cluster_no=3,
        go_interface_instance=None,
        param_set=ref_param_set):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param blank_model_size: the number of uniprots in the blanc model
    :param zoom_range_selector: tuple representing the coverage range for which we would want
     to see the histogram of current distributions
    :param p_val: desired p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used with
     the number of rounds equal to the number
    :param cluster_no: specifies the number of cluster_no we want to have
    :param go_interface_instance:
    :param param_set:
    :return: None if no significant nodes, the node and group characterisitc dictionaries
     otherwise
    """
    if go_interface_instance is None:
        go_interface_instance = get_go_interface_instance(param_set)

    md5_hash = go_interface_instance.md5_hash()

    curr_inf_conf_general = []
    count = 0
    mean_correlation_accumulator = []
    eigenvalues_accumulator = []

    log.info("requested: size: %s, sys_hash: %s, sparse_rounds: %s,",
             blank_model_size, md5_hash,  sparse_rounds)

    log.info("samples found to test against: \t %s",
             annotome_rand_samp.find({'size': blank_model_size,
                                      'sys_hash': md5_hash,
                                      'sparse_rounds': sparse_rounds}).count())

    background_sample = annotome_rand_samp.find(
            {'size': blank_model_size, 'sys_hash': md5_hash, 'sparse_rounds': sparse_rounds})

    for i, sample in enumerate(background_sample):

        if sparse_rounds:
            log.warning('Blank done on sparse rounds. Clustering will not be performed')

        _, node_currents = pickle.loads(sample['currents'])
        tensions = pickle.loads(sample['voltages'])

        if not sparse_rounds:
            _, _, mean_correlations, eigenvalues = perform_clustering(
                tensions, cluster_no, show=False)
        else:
            mean_correlations = np.array([[(0, ), 0, 0]]*cluster_no)
            eigenvalues = np.array([-1]*cluster_no)

        mean_correlation_accumulator.append(np.array(mean_correlations))
        eigenvalues_accumulator.append(eigenvalues)
        dict_system = go_interface_instance.format_node_props(node_currents)
        curr_inf_conf = list(dict_system.values())
        curr_inf_conf_general.append(np.array(curr_inf_conf).T)

        count = i
        if dict_system == {}:
            del mean_correlation_accumulator[-1]
            del curr_inf_conf_general[-1]
            del eigenvalues_accumulator[-1]
            log.critical("exceptional state: nothing in dict_system.")

    # This part declares the pre-operators required for the verification of a
    # real sample

    final = np.concatenate(tuple(curr_inf_conf_general), axis=1)
    final_mean_correlations = np.concatenate(tuple(mean_correlation_accumulator), axis=0).T
    final_eigenvalues = np.concatenate(tuple(eigenvalues_accumulator), axis=0).T

    node_currents = go_interface_instance.node_current
    dict_system = go_interface_instance.format_node_props(node_currents)
    curr_inf_conf_tot = np.array(
        [[int(key)] + list(val) for key, val in dict_system.items()]).T
    go_node_ids, curr_inf_conf = (
        curr_inf_conf_tot[0, :], curr_inf_conf_tot[ (1, 2, 3), :])  # fails here when no significant terms to print
    log.info('blank comparison: %s', curr_inf_conf.shape)
    if not sparse_rounds:
        group2avg_off_diag, _, mean_correlations, eigenvalue = perform_clustering(
            go_interface_instance.UP2UP_voltages, cluster_no, 'GO terms clustering')
    else:
        group2avg_off_diag = np.array([[(0, ), 0, 0]]*cluster_no)
        mean_correlations = np.array([[0, 0]]*cluster_no)
        eigenvalue = np.array([-1]*cluster_no)

    log.info('stats on %s samples', count)

    # TODO: We could and should separate the visualisation from the gaussian
    # estimators computation
    r_nodes, r_groups = show_correlations(
        final, final_mean_correlations, final_eigenvalues,
        zoom_range_selector, curr_inf_conf, mean_correlations.T, eigenvalue.T, count,
        sparse=sparse_rounds, go_interface_instance=go_interface_instance)

    group_char = namedtuple(
        'Group_Char', [
            'UPs', 'num_UPs', 'average_connection', 'p_value'])

    if r_nodes is not None:
        not_random_nodes = [GO_id for GO_id in go_node_ids[r_nodes < p_val].tolist()]

        if not sparse_rounds:
            not_random_groups = np.concatenate(
                (group2avg_off_diag,
                 np.reshape(r_groups, (3, 1))),
                axis=1)[r_groups < p_val].tolist()
            not_random_groups = [group_char(*nr_group)
                                 for nr_group in not_random_groups]
        else:
            not_random_groups = []

        log.debug('not random nodes: %s', not_random_nodes)
        log.debug('bulbs_id_disp_name: %s',
                  list(go_interface_instance.GO2UP_Reachable_nodes.items())[:10])
        # basically the second element below are the nodes that contribute to the information
        #  flow through the node that is considered as non-random

        node_char_list = [
            [int(GO_id), go_interface_instance.GO_Names[GO_id]] +
            dict_system[GO_id] + r_nodes[go_node_ids == float(GO_id)].tolist() +
            [[go_interface_instance.interactome_interface_instance.
                neo4j_id_2_display_name[up_bulbs_id]
             for up_bulbs_id in list(set(go_interface_instance.GO2UP_Reachable_nodes[GO_id]).
                intersection(set(go_interface_instance.analytic_uniprots)))]]
            for GO_id in not_random_nodes]

        return sorted(node_char_list, key=lambda x: x[5]), not_random_groups

    return None, None


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


def auto_analyze(source=None, go_interface_instance=None, processors=3, desired_depth=24,
                 skip_sampling=False, param_set=ref_param_set, output_destination_prefix=''):
    """
    Automatically analyzes the GO annotation of the RNA_seq results.

    :param source:
    :param go_interface_instance:
    :param processors:
    :param desired_depth:
    :param skip_sampling: uses existing mongoDB content without spawning a sampler
    :param param_set:
    """
    if source is None:
        dumplist = undump_object(Dumps.RNA_seq_counts_compare)
    else:
        dumplist = source

    if desired_depth % processors != 0:
        desired_depth = desired_depth // processors + 1
    else:
        desired_depth = desired_depth // processors

    # noinspection PyTypeChecker
    for my_list in dumplist:
        if go_interface_instance is None:
            go_interface_instance = get_go_interface_instance(param_set)

        go_interface_instance.set_uniprot_source(my_list)

        if not skip_sampling:
            log.info("spawning a sampler for %s proteins @ %s compops/sec",
                     len(go_interface_instance.analytic_uniprots), estimated_comp_ops)

        # TODO: restructure to spawn a sampler pool that does not share an object in the Threading
        if len(go_interface_instance.analytic_uniprots) < 200:

            if not skip_sampling:

                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(go_interface_instance.analytic_uniprots),
                         'full',
                         len(go_interface_instance.analytic_uniprots)**2 / estimated_comp_ops /
                         60)

                spawn_sampler_pool(processors,
                                   [len(go_interface_instance.analytic_uniprots)],
                                   [desired_depth],
                                   go_interface_instance=None,
                                   param_set=param_set)

            go_interface_instance.build_extended_conduction_system()
            nr_nodes, nr_groups = compare_to_blank(
                len(go_interface_instance.analytic_uniprots),
                [1100, 1300],
                p_val=0.9,
                go_interface_instance=go_interface_instance,
                param_set=param_set)

        else:
            sampling_depth = max(200 ** 2 // len(go_interface_instance.analytic_uniprots), 5)

            if not skip_sampling:

                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(go_interface_instance.analytic_uniprots),
                         sampling_depth,
                         len(go_interface_instance.analytic_uniprots) * sampling_depth / 2 / 60 / estimated_comp_ops)

                spawn_sampler_pool(processors,
                                   [len(go_interface_instance.analytic_uniprots)],
                                   [desired_depth],
                                   sparse_rounds=sampling_depth,
                                   go_interface_instance=None,
                                   param_set=param_set)

            go_interface_instance.build_extended_conduction_system(sparse_samples=sampling_depth)
            # go_interface_instance.export_conduction_system()
            nr_nodes, nr_groups = compare_to_blank(
                len(go_interface_instance.analytic_uniprots),
                [1100, 1300],
                p_val=0.9, sparse_rounds=sampling_depth,
                go_interface_instance=go_interface_instance,
                param_set=param_set)

        if len(output_destination_prefix) > 0:
            corrected_knowledge_GDF_output = os.path.join(
                os.path.join(Outputs.prefix, output_destination_prefix),
                'GO_Analysis_output.gdf')
            corrected_knowledge_tables_output = os.path.join(
                os.path.join(Outputs.prefix, output_destination_prefix),
                'knowledge_stats.tsv')

        else:
            corrected_knowledge_GDF_output = Outputs.GO_GDF_output
            corrected_knowledge_tables_output = Outputs.knowledge_network_output

        go_interface_instance.export_conduction_system(output_location=corrected_knowledge_GDF_output)

        for group in nr_groups:
            log.info(group)
        log.info('\t NodeID \t Name \t current \t informativity \t confusion_potential \t p_val \t '
                 'UP_list')
        for node in nr_nodes:
            log.info('\t %s \t %s \t %s \t %s \t %s \t %s \t %s', *node)

        with open(corrected_knowledge_tables_output, 'wt') as output:
            # TODO: add an override to the directory of export
            writer = csv_writer(output, delimiter='\t')
            writer.writerow(['NodeID', 'Name', 'current', 'informativity', 'confusion_potential',
                             'p_val', 'UP_list'])
            for node in nr_nodes:
                writer.writerow(node)


if __name__ == "__main__":
    auto_analyze([get_source_bulbs_ids()], processors=3, desired_depth=6, param_set=ref_param_set)
