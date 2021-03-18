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
from typing import Union, List
from collections import defaultdict
from tabulate import tabulate

from bioflow.annotation_network.BioKnowledgeInterface import GeneOntologyInterface
from bioflow.configs.main_configs import estimated_comp_ops, NewOutputs, sparse_analysis_threshold, \
    implicitely_threaded, p_val_cutoff, min_nodes_for_p_val
from bioflow.sample_storage.mongodb import find_annotome_rand_samp, count_annotome_rand_samp
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.io_routines import get_source_bulbs_ids
from bioflow.utils.log_behavior import get_logger
from bioflow.algorithms_bank.flow_significance_evaluation import get_neighboring_degrees, get_p_val_by_gumbel


log = get_logger(__name__)


# CURRENTPASS: [BKI normalization] Those need to go into the user_configs environment
_filter = ['biological_process']
_correlation_factors = (1, 1)
ref_param_set = tuple([_filter, [], (1, 1), True, 3])


# TODO: move ref_param_set to configs
def get_go_interface_instance(param_set=ref_param_set) -> GeneOntologyInterface:
    """
    Generates a Matrix_Knowledge_DB interface for the use in the spawner. If

    :return: a GO_interface object
    """
    go_interface_instance = GeneOntologyInterface(*param_set)
    log.info('Annotation interface created with params %s, %d, %s, %s, %d' %
             (param_set[0], len(param_set[1]), param_set[2], param_set[3], param_set[4]))
    go_interface_instance.fast_load()
    log.info(go_interface_instance.pretty_time())
    return go_interface_instance


# REFACTOR: [SANITY] args puck need to be a named tuple, preferably a typed one
def spawn_sampler(args_puck):
    """
    Spawns a sampler initialized from the default GO_Interface

    :param args_puck: combined list of sample sizes
     and iterations (required for Pool.map usage)
    """
    # log.info('Pool process %d started' % args_puck[-1])

    go_interface_instance = args_puck[3]
    param_set = args_puck[4]

    if go_interface_instance is None:
        go_interface_instance = get_go_interface_instance(param_set)

    sample_size_list = args_puck[0]
    iteration_list = args_puck[1]
    sparse_rounds = args_puck[2]
    pool_no = args_puck[-1]

    go_interface_instance.randomly_sample(
        sample_size_list,
        iteration_list,
        sparse_rounds,
        pool_no=pool_no)

    # log.info('Pool process %d finished' % pool_no)


def spawn_sampler_pool(
        pool_size,
        sample_size_list,
        iterations_list_per_pool,
        sparse_rounds=False,
        go_interface_instance=None,
        param_set=ref_param_set):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param iterations_list_per_pool: number of iterations performing the pooling of the samples
     in each list
    :param sparse_rounds:
    :param go_interface_instance:
    :param param_set: set of parameters configuring the knowledge interface object
    """

    payload = [
        (sample_size_list,
         iterations_list_per_pool,
         sparse_rounds,
         go_interface_instance,
         param_set)]

    payload_list = payload * pool_size
    payload_list = [list(item)+[i] for i, item in enumerate(payload_list)]  # prepare the payload

    global implicitely_threaded

    log.info('Spawning sampler for %s %s' % (payload[0][0], payload[0][1]))
    if not implicitely_threaded:  # TODO: this can be extracted as a shared routine to utils module
        with Pool(processes=pool_size) as pool:  # This is the object we are using to spawn a thread pool
            try:
                log.debug('spawning the sampler with payload %s', payload)
                pool.map(spawn_sampler, payload_list)  # This what we spawn as a sampler
            except Exception as e:
                msg = "{}\n\nOriginal {}".format(e, traceback.format_exc())
                raise type(e)(msg)
            # log.info('Last in-pool flag exiting')
            pool.terminate()  # potential solution to the issue

        log.info('Pool terminated')

    else:
        log.debug('spawning single-thread sampler with payload %s', payload)
        for _payload in payload_list:
            spawn_sampler(_payload)


def samples_scatter_and_hist(background_curr_deg_conf, true_sample_bi_corr_array,
                             save_path: NewOutputs = None, p_values: np.array = None):
    """
    A general function that performs demonstration of an example of random samples of
     the same size as our sample and of our sample and conducts the statistical tests
     on wherther any of nodes or functional groups in our sample are non-random

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
            _filter = p_values < p_val_cutoff
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
    plt.savefig(save_path.knowledge_network_scatterplot)
    plt.clf()


def compare_to_blank(
        blank_model_size: int,
        go_interface_instance: GeneOntologyInterface,
        p_value_cutoff: float = 0.05,
        sparse_rounds: bool = False,
        output_destination: NewOutputs = None):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system

    :param blank_model_size: the number of uniprots in the blanc model
    :param p_value_cutoff: desired p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used with
     the number of rounds equal to the number
    :param cluster_no: specifies the number of cluster_no we want to have
    :param go_interface_instance:
    :param param_set:
    :return: None if no significant nodes, the node and group characterisitc dictionaries
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

    background_sub_array_list = []
    max_sub_array_list = []
    count = 0

    log.info("looking to test against:"
             "\t size: %s \t sys_hash: %s \t sparse_rounds: %s" %
             (blank_model_size, md5_hash,  sparse_rounds))

    log.info("samples found to test against:\t %s" %
             count_annotome_rand_samp({'size': blank_model_size,
                                       'sys_hash': md5_hash,
                                       'sparse_rounds': sparse_rounds}))

    background_sample = find_annotome_rand_samp({'size': blank_model_size,
                                                 'sys_hash': md5_hash,
                                                 'sparse_rounds': sparse_rounds})

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
        if dict_system == {}:
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

    curr_inf_conf_tot = np.array([[int(key)] + list(val) for key, val in dict_system.items()]).T

    go_node_ids, query_array = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1, 2, 3), :])

    # new p-values computation
    background_density = kde_compute(background_array[(1, 0), :], 50, count)
    base_bi_corr = background_array[(0, 1), :]

    r_rels = []
    r_std_nodes = []

    degrees = np.unique(query_array[2, :])

    combined_p_vals = np.ones_like(query_array[2, :])

    for degree in degrees.tolist():  # TODO: there is currently a logic where the pb fails if
        # there is elements with edge number in the analysis set that were not found in the test
        # sets (which is not unexpected for GO terms)
        _filter = query_array[2, :] == degree

        entry = query_array[:, _filter]
        background_set = background_array[:, background_array[2, :] == degree]

        max_current_per_run = get_neighboring_degrees(degree,
                                                      max_array,
                                                      min_nodes=min_nodes_for_p_val)

        p_vals = get_p_val_by_gumbel(entry, max_current_per_run)
        combined_p_vals[_filter] = p_vals

    samples_scatter_and_hist(background_array, query_array,
                             save_path=output_destination,
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
              list(go_interface_instance.GO2UP_Reachable_nodes.items())[:10])

    node_char_list = [
            [int(GO_id),
             go_interface_instance.GO_Names[GO_id]] +
            dict_system[GO_id] +
            r_nodes[go_node_ids == float(GO_id)].tolist() +
            [[go_interface_instance.UP_Names[up_bulbs_id][1]
              for up_bulbs_id in list(set(go_interface_instance.GO2UP_Reachable_nodes[GO_id]).
                                      intersection(set(go_interface_instance.active_up_sample)))]]
            for GO_id in not_random_nodes]

    nodes_dict = np.hstack((go_node_ids[:, np.newaxis],
                            r_nodes[:, np.newaxis],
                            r_rels[:, np.newaxis],
                            r_std_nodes[:, np.newaxis]))
    nodes_dict = dict((node[0], (node[1], node[2], node[3])) for node in nodes_dict.tolist())
    nodes_dict = defaultdict(lambda: (1., 0., 0.), nodes_dict)  # corresponds to the cases of super low flow - never significant

    # TODO: pull the groups corresponding to non-random associations.
    # => Will not implement, it's already done by Gephi

    return sorted(node_char_list, key=lambda x: x[5]), nodes_dict


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


# TODO: [weighted inputs] add support for a dict as source_list, not only list
def auto_analyze(source_list,
                 output_destinations_list: Union[List[str], None] = None,
                 go_interface_instance=None,
                 processors=3,
                 desired_depth=24,
                 skip_sampling=False,
                 param_set=ref_param_set,
                 p_value_cutoff: float = -1,
                 ) -> None:
    """
    Automatically analyzes the GO annotation of the RNA_seq results.

    :param source_list:
    :param go_interface_instance:
    :param processors:
    :param desired_depth:
    :param skip_sampling: uses existing mongoDB content without spawning a sampler
    :param param_set:
    """

    # Multiple re-spawns of threaded processing are incompatbile with scikits.sparse.cholmod
    if len(source_list) > 1:
        global implicitely_threaded
        implicitely_threaded = True

    if len(output_destinations_list) != len(source_list):
        log.warning('Output destination list has %d elements, whereas %d sources were supplied. '
                    'Falling back to default output structure')
        output_destinations_list = None

    if output_destinations_list is None:
        output_destinations_list = list(range(len(source_list)))

    if processors == 0:
        processors = psutil.cpu_count() - 1
        log.info("Setting processor count to default: %s" % processors)

    if desired_depth % processors != 0:
        desired_depth = desired_depth // processors + 1
    else:
        desired_depth = desired_depth // processors

    if p_value_cutoff < 0:
        p_value_cutoff = p_val_cutoff

    for hits_list, output_destination in zip(source_list, output_destinations_list):

        log.info('Auto analyzing list of interest: %s', len(hits_list))

        outputs_subdirs = NewOutputs(output_destination)

        if go_interface_instance is None:
            go_interface_instance = get_go_interface_instance(param_set)

        go_interface_instance.set_uniprot_source(hits_list)

        # go_intraface instance sets the background if it has one by itself.
        # the background setting, however, works in a different way => check that.

        if not skip_sampling:
            log.info("spawning a sampler for %s proteins @ %s compops/sec",
                     len(go_interface_instance.active_up_sample), estimated_comp_ops)

        if len(go_interface_instance.active_up_sample) < sparse_analysis_threshold:

            if not skip_sampling:

                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(go_interface_instance.active_up_sample),
                         'full',
                         len(go_interface_instance.active_up_sample) ** 2 / estimated_comp_ops /
                         60)

                spawn_sampler_pool(processors,
                                   [len(go_interface_instance.active_up_sample)],
                                   [desired_depth],
                                   go_interface_instance=None,
                                   param_set=param_set)

            go_interface_instance.compute_current_and_potentials()

            nr_nodes, p_val_dict = compare_to_blank(
                len(go_interface_instance.active_up_sample),
                go_interface_instance,
                p_value_cutoff=p_value_cutoff,
                output_destination=outputs_subdirs)

        # sparse analysis
        else:
            ceiling = min(205, len(go_interface_instance.active_up_sample))
            sampling_depth = max((ceiling - 5) ** 2 // len(go_interface_instance.active_up_sample), 5)

            if not skip_sampling:

                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(go_interface_instance.active_up_sample),
                         sampling_depth,
                         len(go_interface_instance.active_up_sample) *
                         sampling_depth / 2 / 60 / estimated_comp_ops)

                spawn_sampler_pool(processors,
                                   [len(go_interface_instance.active_up_sample)],
                                   [desired_depth],
                                   sparse_rounds=sampling_depth,
                                   go_interface_instance=None,
                                   param_set=param_set)

            go_interface_instance.compute_current_and_potentials(sparse_samples=sampling_depth)

            # go_interface_instance.export_conduction_system()
            nr_nodes, p_val_dict = compare_to_blank(
                len(go_interface_instance.active_up_sample),
                go_interface_instance,
                p_value_cutoff=p_value_cutoff,
                sparse_rounds=sampling_depth,
                output_destination=outputs_subdirs)

        go_interface_instance.export_conduction_system(p_val_dict,
                                                       output_location=outputs_subdirs.GO_GDF_output)


        # # old results print-out
        # log.info('\t NodeID \t Name \t current \t informativity \t confusion_potential \t p_val \t '
        #          'UP_list')
        #
        # for node in nr_nodes:
        #     log.info('\t %s \t %s \t %.3g \t %.3g \t %d \t %.3g \t %s', *node)

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


if __name__ == "__main__":
    auto_analyze([get_source_bulbs_ids()], processors=3, desired_depth=6, param_set=ref_param_set)
