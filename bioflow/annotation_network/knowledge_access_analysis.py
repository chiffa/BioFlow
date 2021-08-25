"""
Set of methods responsible for knowledge analysis
"""
import pickle
from multiprocessing import Pool
import traceback
import psutil
import numpy as np
from matplotlib import pyplot as plt
from csv import writer as csv_writer
from collections import defaultdict
from tabulate import tabulate

from typing import Any, Union, TypeVar, NewType, Tuple, List

from bioflow.annotation_network.BioKnowledgeInterface import GeneOntologyInterface
from bioflow.configs.main_configs import estimated_comp_ops, NewOutputs, sparse_analysis_threshold, \
    implicitely_threaded, default_p_val_cutoff, min_nodes_for_p_val, default_background_samples
from bioflow.sample_storage.mongodb import find_annotome_rand_samp, count_annotome_rand_samp
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.io_routines import get_source_bulbs_ids
from bioflow.utils.log_behavior import get_logger
from bioflow.algorithms_bank.flow_significance_evaluation import get_neighboring_degrees, \
    get_p_val_by_gumbel
from bioflow.algorithms_bank.clustering_routines import compute_tension_clustering
import bioflow.algorithms_bank.sampling_policies as sampling_policies

log = get_logger(__name__)


def get_go_interface_instance(background: Union[List[int], List[Tuple[int, float]]] = ()) -> \
        GeneOntologyInterface:
    """
    Generates a Matrix_Knowledge_DB interface for the use in the spawner. If

    :return: a GO_interface object
    """
    go_interface_instance = GeneOntologyInterface(background=background)
    go_interface_instance.fast_load()
    log.info('Annotation interface created with background length %d' %
             len(go_interface_instance._background))
    return go_interface_instance


# REFACTOR: [MAINTAINABILITY] args puck need to be a named tuple, preferably a typed one
def _spawn_sampler(args_puck):
    """
    Spawns a sampler initialized from the default GO_Interface

    :param args_puck: combined list of sample sets to match, expected random samples and other
        parameters (packing required to start the pool)
    """
    # (sample_sets_to_match,        0
    # sample_depth,                 1
    # sparse_rounds,                2
    # background_set,               3
    # forced_go_interface,          4
    # sampling_policy,              5
    # sampling_options)             6
    # pool_no (added by map)        7

    background_set_arg = args_puck[3]

    if args_puck[4] is not None:
        go_interface_instance = args_puck[4]
    else:
        go_interface_instance = get_go_interface_instance(background_set_arg)

    hits_list, sec_list = args_puck[0]
    go_interface_instance.set_flow_sources(hits_list, sec_list)

    iterations = args_puck[1]
    sparse_rounds = args_puck[2]

    sampling_policy = args_puck[5]
    sampling_options = args_puck[6]

    pool_no = args_puck[-1]

    go_interface_instance.reset_thread_hex()
    go_interface_instance.randomly_sample(
        iterations,
        sparse_rounds=sparse_rounds,
        pool_no=pool_no,
        sampling_policy=sampling_policy,
        optional_sampling_param=sampling_options
    )


def _spawn_sampler_pool(
        pool_size,
        sample_sets_to_match,
        sample_depth,
        sparse_rounds,
        background_set,
        forced_go_interface,
        sampling_policy,
        sampling_options):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_sets_to_match: size of the sample list
    :param sample_depth: number of random samples we look to generate performing the pooling of the
        samples in each list
    :param sparse_rounds: number of sparse rounds to run (or False if sampling is dense)
    :param background_set: set of node ids that are to be sampled from
    :param forced_go_interface: a provided BioKnowledgeInterface that contains sets to imitate
    :param sampling_policy: sampling policy to be employed
    :param sampling_options: options to the sampling policy
    """
    global implicitely_threaded

    if not implicitely_threaded:
        samples_remaining = sample_depth % pool_size

        if samples_remaining > 0:
            sample_depth = sample_depth // pool_size + 1
        else:
            sample_depth = sample_depth // pool_size

    payload = [(sample_sets_to_match,
                sample_depth,
                sparse_rounds,
                background_set,
                forced_go_interface,
                sampling_policy,
                sampling_options)]

    payload_list = payload * pool_size
    payload_list = [list(item)+[i] for i, item in enumerate(payload_list)]  # prepare the payload

    # log.info('Spawning sampler for %s %s' % (payload[0][0], payload[0][1]))
    if not implicitely_threaded:
        with Pool(processes=pool_size) as pool:  # This is the object we are using to spawn a thread pool
            try:
                log.debug('spawning the sampler with payload %s', payload)
                pool.map(_spawn_sampler, payload_list)  # This what we spawn as a sampler
                # KNOWNBUG: hangs with no message upon a second start attempt in Interactome
                #  analysis due to cholmod
            except Exception as e:
                msg = "{}\n\nOriginal {}".format(e, traceback.format_exc())
                raise type(e)(msg)
            # log.info('Last in-pool flag exiting')
            pool.terminate()  # potential solution to the issue

        # log.info('Pool terminated')

    else:
        log.debug('spawning single-thread sampler with payload %s', payload)
        for _payload in payload_list:
            _spawn_sampler(_payload)


def samples_scatter_and_hist(background_curr_deg_conf, true_sample_bi_corr_array,
                             save_path: str = None, p_values: np.array = None):
    """
    A general function that performs demonstration of an example of random samples of
     the same size as our sample and of our sample and conducts the statistical tests
     on whether any of nodes or functional groups in our sample are non-random

    :param background_curr_deg_conf: [[current, informativity, confusion_potential], ...] -
    characteristics of the random samples
    :param true_sample_bi_corr_array: [[current, informativity, confusion_potential], ...] -
    characteristics of the true sample. If none, nothing happens
    :param save_path: where the thing will be saved
    :param p_values: p-value map that will be used to save things after the analysis
    :return: None
    """

    fig = plt.figure()
    fig.set_size_inches(30, 20)

    # bivect: [0, :] - current; [1, :] - informativity

    plt.subplot(211)
    plt.title('current through nodes')

    bins = np.linspace(
        background_curr_deg_conf[0, :].min(),
        background_curr_deg_conf[0, :].max(), 100)

    if true_sample_bi_corr_array is not None:
        bins = np.linspace(min(background_curr_deg_conf[0, :].min(),
                               true_sample_bi_corr_array[0, :].min()),
                           max(background_curr_deg_conf[0, :].max(),
                               true_sample_bi_corr_array[0, :].max()),
                           100)

    plt.hist(background_curr_deg_conf[0, :],
             bins=bins, histtype='step', log=True, color='b')

    if true_sample_bi_corr_array is not None:
        plt.hist(true_sample_bi_corr_array[0, :],
                 bins=bins, histtype='step', log=True, color='r')


    plt.subplot(212)
    plt.scatter(background_curr_deg_conf[1, :],
                background_curr_deg_conf[0, :], color='b', alpha=0.1)

    if true_sample_bi_corr_array is not None:
        if p_values is not None:
            _filter = p_values < default_p_val_cutoff
            anti_filter = np.logical_not(_filter)
            plt.scatter(true_sample_bi_corr_array[1, anti_filter],
                        true_sample_bi_corr_array[0, anti_filter],
                        color='gray', alpha=0.25)

            plt.scatter(true_sample_bi_corr_array[1, _filter],
                        true_sample_bi_corr_array[0, _filter],
                        color='r', alpha=0.7)

        else:
            plt.scatter(true_sample_bi_corr_array[1, :],
                        true_sample_bi_corr_array[0, :],
                        color='r', alpha=0.5)

    # plt.show()
    if save_path is not None:
        plt.savefig(save_path)

    plt.clf()


def clustering_analysis_complement(go_interface_instance: GeneOntologyInterface,
                                   p_val_cutoff: float = 0.05,
                                   sparse_rounds: int = -1,
                                   output_destination: NewOutputs = None,
                                   random_sampling_method=sampling_policies.matched_sampling,
                                   random_sampling_option='exact',
                                   ) -> list:

    def get_max_for_each_degree(min_inf_flow, clust_size):
        # print('debug max_array_shape:', str(sample_sub_arrray.shape))
        degrees = np.unique(clust_size)
        _max_array = []

        for degree in degrees:
            l_filter = clust_size == degree
            _max_array.append([min_inf_flow[l_filter].max(), degree])

        m_arr = np.array(_max_array)
        return m_arr.T

    if sparse_rounds > 0:
        return []

    if go_interface_instance is None or go_interface_instance.node_current == {}:
        raise Exception("tried to perform clustering complement analysis on an empty interface "
                        "instance")

    md5_hash = go_interface_instance.md5_hash()
    active_sample_hash = go_interface_instance.active_sample_md5_hash(sparse_rounds)

    background_sub_array_list = []
    max_sub_array_list = []

    count = 0

    log.info("GO clustering complement looking to test against:\n"
             "\t target_hash: %s \t sys_hash: %s \n"
             "random sampled according to %s/%s" %
             (active_sample_hash, md5_hash, random_sampling_method.__name__, random_sampling_option))

    samples_to_test_against = count_annotome_rand_samp({
                                          'active_sample_hash': active_sample_hash,
                                          'sys_hash': md5_hash,
                                          'sampling_policy': random_sampling_method.__name__,
                                          'sampling_policy_options': random_sampling_option})

    if samples_to_test_against == 0:
        raise Exception('No samples found to test against. '
                        'There is likely a discrepancy between the parameters used for sampling '
                        'and parameters used for the life sample analysis')

    log.info("samples found to test against:\t %d" % samples_to_test_against)

    background_samples = find_annotome_rand_samp({
                                          'active_sample_hash': active_sample_hash,
                                          'sys_hash': md5_hash,
                                          'sampling_policy': random_sampling_method.__name__,
                                          'sampling_policy_options': random_sampling_option})

    for i, sample in enumerate(background_samples):

        voltages = pickle.loads(sample['voltages'])

        _, min_clust_inf_flow, clust_size = compute_tension_clustering(voltages, random_sample=True)

        back_arr = np.vstack([min_clust_inf_flow, clust_size])
        background_sub_array_list.append(back_arr)

        max_arr = get_max_for_each_degree(min_clust_inf_flow, clust_size)
        max_sub_array_list.append(max_arr)
        count = i

    background_array = np.concatenate(tuple(background_sub_array_list), axis=1)
    max_array = np.concatenate(tuple(max_sub_array_list), axis=1)

    voltages = go_interface_instance.UP2UP_voltages
    clusters, min_clust_inf_flow, clust_size = compute_tension_clustering(voltages,
                                                                          random_sample=False)

    query_array = np.vstack([min_clust_inf_flow, clust_size])

    clust_sizes = np.unique(clust_size)

    combined_p_vals = np.ones_like(min_clust_inf_flow)

    for cluster_size in clust_sizes.tolist():
        _filter = clust_size == cluster_size
        entry = min_clust_inf_flow[_filter]

        min_clust_inf_flow_per_run = get_neighboring_degrees(cluster_size,
                                                             max_array,
                                                             min_nodes=min_nodes_for_p_val)

        p_vals = get_p_val_by_gumbel(entry, min_clust_inf_flow_per_run)
        combined_p_vals[_filter] = p_vals

    samples_scatter_and_hist(background_array, query_array,
                             save_path=output_destination.knowledge_clusters_scatterplot,
                             # to save
                             p_values=combined_p_vals)

    not_random_clusters = [(cluster, i)
                           for i, cluster
                           in enumerate(clusters)
                           if combined_p_vals[i] < p_val_cutoff]

    # rough render of the results

    cluster_entries = []

    for cluster, cluster_no in not_random_clusters:
        cluster_entries.append([cluster_no,
                                combined_p_vals[cluster_no],
                                clust_size[cluster_no],
                                min_clust_inf_flow[cluster_no],
                                []])
        for internal_id in cluster:
            leg_id, f_name = go_interface_instance.up_neo4j_id_2_leg_id_disp_name[internal_id]
            legacy_id = leg_id
            node_type = 'NA (UP)'
            full_name = f_name
            cluster_entries[-1][-1].append([internal_id, legacy_id, node_type, full_name])

    # sort by p-value:
    cluster_entries = sorted(cluster_entries, key=lambda x: x[1]    )

    return cluster_entries


def compare_to_blank(
        go_interface_instance: GeneOntologyInterface,
        p_value_cutoff: float = 0.05,
        sparse_rounds: int = -1,
        output_destination: NewOutputs = None,
        random_sampling_method=sampling_policies.matched_sampling,
        random_sampling_option='exact',
        ) -> Tuple[list, dict]:
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param go_interface_instance:
    :param p_value_cutoff: desired p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used with
        the number of rounds equal to the number
    :param output_destination: configs object from main_configs, specifying where the results
        will be saved
    :param random_sampling_method: sampling policy used
    :param random_sampling_option: sampling policy optional argument
    :return: None if no significant nodes, the node and group characteristic dictionaries
     otherwise
    """

    def get_max_for_each_degree(sample_sub_arrray):
        # this will work only if we are using confusion potential (which is the # of
        #  nodes a term annotates)
        # print('debug max_array_shape:', str(sample_sub_arrray.shape))
        degrees = np.unique(sample_sub_arrray[2, :])
        max_array = []

        for degree in degrees:
            filter = sample_sub_arrray[2, :] == degree
            max_array.append([sample_sub_arrray[0, filter].max(), degree])

        m_arr = np.array(max_array)
        return m_arr.T

    if go_interface_instance is None or go_interface_instance.node_current == {}:
        raise Exception("tried to compare to blanc an empty interface instance")

    md5_hash = go_interface_instance.md5_hash()
    active_sample_hash = go_interface_instance.active_sample_md5_hash(sparse_rounds)

    background_sub_array_list = []
    max_sub_array_list = []
    count = 0

    log.info("looking to test against:"
             "\t target_hash: %s \t sys_hash: %s \n"
             "random sampled according to %s/%s" %
             (active_sample_hash, md5_hash, random_sampling_method.__name__, random_sampling_option))

    samples_to_test_against = count_annotome_rand_samp({
                                          'active_sample_hash': active_sample_hash,
                                          'sys_hash': md5_hash,
                                          'sampling_policy': random_sampling_method.__name__,
                                          'sampling_policy_options': random_sampling_option})

    if samples_to_test_against == 0:
        raise Exception('No samples found to test against. '
                        'There is likely a discrepancy between the parameters used for sampling '
                        'and parameters used for the life sample analysis')

    log.info("samples found to test against:\t %d" % samples_to_test_against)

    background_sample = find_annotome_rand_samp({
                                          'active_sample_hash': active_sample_hash,
                                          'sys_hash': md5_hash,
                                          'sampling_policy': random_sampling_method.__name__,
                                          'sampling_policy_options': random_sampling_option})

    for i, sample in enumerate(background_sample):

        _, node_currents = pickle.loads(sample['currents'])

        dict_system = go_interface_instance.format_node_props(node_currents)
        background_sub_array = list(dict_system.values())

        if np.array(background_sub_array).T.shape[0] < 2:
            log.info(background_sub_array)
            continue

        background_sub_array_list.append(np.array(background_sub_array).T)

        max_arr = get_max_for_each_degree(np.array(background_sub_array).T)
        max_sub_array_list.append(max_arr)

        count = i

        if dict_system == {}:  # ?????
            del background_sub_array_list[-1]
            del max_sub_array_list[-1]
            log.critical("exceptional state: nothing in dict_system. Attempting to ignore")

    # This part declares the pre-operators required for the verification of a
    # real sample

    background_array = np.concatenate(tuple(background_sub_array_list), axis=1)
    max_array = np.concatenate(tuple(max_sub_array_list), axis=1)

    # final = np.concatenate(tuple(background_sub_array_list), axis=1)
    # final_mean_correlations = np.concatenate(tuple(mean_correlation_accumulator), axis=0).T
    # final_eigenvalues = np.concatenate(tuple(eigenvalues_accumulator), axis=0).T

    node_currents = go_interface_instance.node_current
    dict_system = go_interface_instance.format_node_props(node_currents)

    curr_inf_conf_tot = np.array([[int(key)] + list(val) for key, val in list(dict_system.items())]).T

    go_node_ids, query_array = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1, 2, 3), :])

    log.info("stats on  %s samples" % count)
    # new p-values computation
    background_density = kde_compute(background_array[(1, 0), :], 50, count)
    base_bi_corr = background_array[(0, 1), :]

    r_rels = []
    r_std_nodes = []

    degrees = np.unique(query_array[2, :])

    combined_p_vals = np.ones_like(query_array[2, :])

    for degree in degrees.tolist():
        _filter = query_array[2, :] == degree

        entry = query_array[0, _filter]
        background_set = background_array[:, background_array[2, :] == degree]

        # REFACTOR: [maintenability] this part is too coupled. we should factor it out
        max_current_per_run = get_neighboring_degrees(degree,
                                                      max_array,
                                                      min_nodes=min_nodes_for_p_val)

        p_vals = get_p_val_by_gumbel(entry, max_current_per_run)
        combined_p_vals[_filter] = p_vals

    samples_scatter_and_hist(background_array, query_array,
                             save_path=output_destination.knowledge_network_scatterplot,
                             p_values=combined_p_vals)

    r_nodes = background_density(query_array[(1, 0), :])  # legacy - unused now
    r_nodes = combined_p_vals

    for point in query_array.T:
        selector = np.logical_and(base_bi_corr[1, :] > point[1]*0.9, base_bi_corr[1, :] < point[1]*1.1)
        r_rels.append(point[0] / np.mean(base_bi_corr[0, selector]))
        r_std_nodes.append((point[0] - np.mean(base_bi_corr[0, selector])) / np.std(base_bi_corr[0,
                                                                                              selector]))

    r_rels = np.array(r_rels)
    r_std_nodes = np.array(r_std_nodes)

    not_random_nodes = [node_id for node_id in go_node_ids[r_nodes < p_value_cutoff].tolist()]

    log.debug('debug, not random nodes: %s', not_random_nodes)
    log.debug('debug bulbs_id_disp_name: %s',
              list(go_interface_instance._limiter_go_2_up_reachable_nodes.items())[:10])

    node_char_list = [
        [int(nr_node_id),
         go_interface_instance.neo4j_id_2_display_name[nr_node_id]] +
        dict_system[nr_node_id] + r_nodes[go_node_ids == float(nr_node_id)].tolist() +
        [[go_interface_instance.up_neo4j_id_2_leg_id_disp_name[up_bulbs_id][1]  # UP name
              for up_bulbs_id
              in list(set(go_interface_instance._limiter_go_2_up_reachable_nodes[nr_node_id]).
                      intersection(set(go_interface_instance._active_up_sample)))]]
        for nr_node_id in not_random_nodes]

    nodes_dict = np.hstack((go_node_ids[:, np.newaxis],
                            r_nodes[:, np.newaxis],
                            r_rels[:, np.newaxis],
                            r_std_nodes[:, np.newaxis]))
    nodes_dict = dict((node[0], (node[1], node[2], node[3])) for node in nodes_dict.tolist())
    nodes_dict = defaultdict(lambda: (1., 0., 0.), nodes_dict)  # corresponds to the cases of super low flow - never significant

    return sorted(node_char_list, key=lambda x: x[5]), nodes_dict


def auto_analyze(source_list: List[Union[List[int], List[Tuple[int, float]]]],
                 secondary_source_list: List[Union[List[int], List[Tuple[int, float]], None]] = None,
                 output_destinations_list: Union[List[str], None] = None,
                 random_samples_to_test_against=default_background_samples,
                 processors: int = 0,
                 background_list: List[Union[List[int], List[Tuple[int, float]], None]] = None,
                 skip_sampling: bool = False,
                 cluster: bool = True,
                 p_value_cutoff: float = -1,
                 sampling_policy=sampling_policies.matched_sampling,
                 sampling_policy_options='exact',
                 explicit_interface=None,
                 ) -> None:
    """
    Automatically analyzes the GO annotation of the experimental hit lists

    :param source_list: python list of hits for each condition
    :param secondary_source_list: secondary list to which calculate the flow from hist, if needed
    :param output_destinations_list: list of names for each condition
    :param random_samples_to_test_against: total samples we would like to compare each set of hits with
    :param processors: number of processes that will be loaded. as a rule of thumb,
        for max performance, use N-1 processors, where N is the number of physical cores on the
        machine, which is the default
    :param background_list:  list of physical entities that an experimental method can retrieve,
        optionally with weights indicating the likelihood of retrieval at random
    :param skip_sampling: if true, will skip background sampling step
    :param cluster: if set to true, will perform the clustering analysis complement
    :param p_value_cutoff: highest p_value up to which to report the results
    :param sampling_policy: sampling policy used
    :param sampling_policy_options: sampling policy optional argument
    :param explicit_interface: an explicit BioKnowledgeInterface instance in case any of the deep
        defaults (eg flow calculation function) are modified
    :return:
    """
    # Multiple re-spawns of threaded processing are incompatible with scikits.sparse.cholmod

    if background_list is None:
        background_list = []

    if len(source_list) > 1:
        global implicitely_threaded
        implicitely_threaded = True

    if output_destinations_list is None:
        output_destinations_list = list(range(len(source_list)))
        output_destinations_list = [str(_item) for _item in output_destinations_list]

    if len(output_destinations_list) != len(source_list):  # we are not calling len on None
        log.warning('Output destination list has %d elements, whereas %d sources were supplied. '
                    'Falling back to default output structure')
        output_destinations_list = list(range(len(source_list)))
        output_destinations_list = [str(_item) for _item in output_destinations_list]

    if processors == 0:
        processors = psutil.cpu_count() - 1
        log.info("Setting processor count to default: %s" % processors)

    if p_value_cutoff <= 0:
        p_value_cutoff = default_p_val_cutoff

    if secondary_source_list is None:
        secondary_source_list = [None] * len(source_list)

    for hits_list, sec_list, output_destination in zip(source_list, secondary_source_list,
                                                       output_destinations_list):

        if hits_list is None or len(hits_list) < 2:
            log.warning('hits list for destination %s contains less than two items: (%s).'
                        'Skipping the analysis' % (output_destination, hits_list))
            continue

        log.debug('debug 2 : %s, %s' % (hits_list, sec_list))

        prim_len, prim_shape, _, sec_len, sec_shape, _, _, _ = \
            sampling_policies.characterize_flow_parameters(hits_list, sec_list, False)

        log.info('Auto analyzing hits list of shapes: %d/%d; %d/%d' %
                 (prim_len, prim_shape, sec_len, sec_shape))

        outputs_subdirs = NewOutputs(output_destination)

        go_interface = get_go_interface_instance(background=background_list)

        if explicit_interface is not None:
            go_interface = explicit_interface

        go_interface.set_flow_sources(hits_list, sec_list)
        total_ops = go_interface.evaluate_ops()
        sparse_rounds = go_interface.reduce_ops(sparse_analysis_threshold**2)

        if sparse_rounds > 0:
            log.info('estimated ops for dense sampling would have been %.1f, '
                     'which is more than the threshold (%d^2). '
                     'Using sparse sampling with %d rounds' % (total_ops,
                                                               sparse_analysis_threshold,
                                                               sparse_rounds))

        else:
            log.info('estimated ops for dense sampling %.1f' % (total_ops))

        md5_hash = go_interface.md5_hash()
        active_sample_hash = go_interface.active_sample_md5_hash(sparse_rounds)

        in_storage = count_annotome_rand_samp({'active_sample_hash': active_sample_hash,
                                               'sys_hash': md5_hash,
                                               'sampling_policy': sampling_policy.__name__,
                                               'sampling_policy_options': sampling_policy_options})

        if in_storage > random_samples_to_test_against:
            log.info("%d suitable random samples found in storage for %d desired. Skipping "
                     "sampling" % (in_storage, random_samples_to_test_against))
            skip_sampling = True

        else:
            log.info("%d suitable random samples found in storage for %d desired. Sampling %d" %
                     (in_storage, random_samples_to_test_against, random_samples_to_test_against - in_storage))
            random_samples_to_test_against = random_samples_to_test_against - in_storage

        if not skip_sampling:
            _spawn_sampler_pool(processors,
                                (hits_list, sec_list),
                                random_samples_to_test_against,
                                sparse_rounds=sparse_rounds,
                                background_set=background_list,
                                forced_go_interface=explicit_interface,
                                sampling_policy=sampling_policy,
                                sampling_options=sampling_policy_options)

        go_interface.compute_current_and_potentials(sparse_rounds=sparse_rounds)

        nr_nodes, p_val_dict = compare_to_blank(
            go_interface,
            p_value_cutoff=p_value_cutoff,
            sparse_rounds=sparse_rounds,
            output_destination=outputs_subdirs,
            random_sampling_method=sampling_policy,
            random_sampling_option=sampling_policy_options
        )

        go_interface.export_conduction_system(p_val_dict,
                                                    output_location=outputs_subdirs.GO_GDF_output)

        with open(outputs_subdirs.knowledge_network_output, 'wt') as output:
            writer = csv_writer(output, delimiter='\t')
            writer.writerow(['NodeID', 'Name', 'current', 'informativity', 'confusion_potential',
                             'p_val', 'UP_list'])
            for node in nr_nodes:
                writer.writerow(node)

        # using tabulate

        headers = ['NodeID', 'Name', 'current', 'informativity', 'confusion_potential', 'p_val',
                   'UP_list']

        print(tabulate(nr_nodes, headers, tablefmt='simple', floatfmt=".3g"))

        if cluster:
            cluster_entries = clustering_analysis_complement(
                go_interface,
                p_val_cutoff=p_value_cutoff,
                sparse_rounds=sparse_rounds,
                output_destination=outputs_subdirs,
                random_sampling_method=sampling_policy,
                random_sampling_option=sampling_policy_options
            )

            with open(outputs_subdirs.knowledge_clusters_output, 'wt') as output:
                writer = csv_writer(output, delimiter='\t')
                for cluster in cluster_entries:
                    writer.writerow(['cluster_no', 'cluster_p_val', 'cluster_size',
                                     'min_cluster_info_flow'])
                    writer.writerow(cluster[:-1])
                    writer.writerow(['', 'id', 'legacy id', 'type', 'display name'])
                    for node in cluster[-1]:
                        writer.writerow([''] + node)

            # using tabulate output to stdout:

            cluster_headers = ['cluster_no', 'cluster_p_val', 'cluster_size', 'min_cluster_info_flow']
            nodes_list_headers = ['', 'id', 'legacy id', 'type', 'display name']

            for cluster in cluster_entries:
                print(tabulate([cluster[:-1]], cluster_headers, tablefmt='simple', floatfmt=".3g"))
                print(tabulate(cluster[-1], nodes_list_headers, tablefmt='simple', floatfmt=".3g"))


if __name__ == "__main__":
    source, sec_source = get_source_bulbs_ids()
    auto_analyze(source, processors=3, random_samples_to_test_against=6, background_list=[])
