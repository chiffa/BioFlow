"""
New analytical routines for the interactome
"""
import pickle
from collections import namedtuple
from csv import reader
from csv import writer as csv_writer
from multiprocessing import Pool
from collections import defaultdict
import traceback
from pprint import pprint
import os
import psutil
from typing import Any, Union, TypeVar, NewType, Tuple, List
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gumbel_r
from tabulate import tabulate

from bioflow.configs.main_configs import Dumps, estimated_comp_ops, NewOutputs, \
    sparse_analysis_threshold, implicitely_threaded, default_p_val_cutoff, min_nodes_for_p_val
from bioflow.sample_storage.mongodb import find_interactome_rand_samp, count_interactome_rand_samp
from bioflow.configs.main_configs import output_location
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.log_behavior import get_logger
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids
from bioflow.utils.general_utils.high_level_os_io import mkdir_recursive
from bioflow.algorithms_bank.flow_significance_evaluation import get_neighboring_degrees,\
    active_default_p_value_policy


log = get_logger(__name__)


def get_interactome_interface(background_up_ids=()) -> InteractomeInterface:
    """
    Retrieves an "InteractomeInterface" object

    :return:
    """
    interactome_interface_instance = InteractomeInterface(background_up_ids=background_up_ids)
    interactome_interface_instance.fast_load()
    log.debug("get_interactome state e_p_u_b_i length: %s",
              len(interactome_interface_instance.active_up_sample))
    log.info("interactome interface loaded in %s" % interactome_interface_instance.pretty_time())
    # is the case now
    return interactome_interface_instance


def spawn_sampler(args_puck):
    """
    Spawns a sampler initialized from the default GO_Interface.

    :param args_puck: combined list of sample sizes, iterations, background sets, and sparse
    sampling argument
    """
    # log.info('Pool process %d started' % args_puck[-1])

    background_set_arg = args_puck[3]
    interactome_interface_instance = get_interactome_interface(background_set_arg)

    sample_size_list = args_puck[0]
    iteration_list = args_puck[1]
    sparse_rounds = args_puck[2]
    pool_no = args_puck[-1]   # TODO: switch over to PID here

    interactome_interface_instance.reset_thread_hex()
    interactome_interface_instance.randomly_sample(
        sample_size_list,
        iteration_list,
        sparse_rounds,
        pool_no=pool_no
    )


def spawn_sampler_pool(
        pool_size,
        sample_size_list,
        interaction_list_per_pool,
        background_set,
        sparse_rounds=False):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param interaction_list_per_pool: number of iterations performing the pooling of the samples
     in each list
    :param sparse_rounds: number of sparse rounds to run (or False if sampling is dense)
    :param background_set: set of node ids that are to be sampled from
    """
    payload = [
            (sample_size_list,
             interaction_list_per_pool,
             sparse_rounds,
             background_set)]

    payload_list = payload * pool_size
    payload_list = [list(item)+[i] for i, item in enumerate(payload_list)]  # prepare the payload

    global implicitely_threaded

    if not implicitely_threaded:
        with Pool(processes=pool_size) as pool:  # This is the object we are using to spawn a thread pool
            try:
                log.debug('spawning the sampler with payload %s', payload)
                pool.map(spawn_sampler, payload_list)  # This what we spawn as a sampler
                # KNOWNBUG: hangs with no message upon a second start attempt in Interactome
                #  analysis due to cholmod
            except Exception as e:
                msg = "{}\n\nOriginal {}".format(e, traceback.format_exc())
                raise type(e)(msg)
            # log.info('Last in-pool flag exiting')
            pool.terminate()

        # log.info('Pool terminated')

    else:
        log.debug('spawning single-thread sampler with payload %s', payload)
        for _payload in payload_list:
            spawn_sampler(_payload)


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
    plt.savefig(save_path.interactome_network_scatterplot)
    plt.clf()


def compare_to_blank(blank_model_size: int,
                     interactome_interface_instance: InteractomeInterface,
                     p_val: float = 0.05,
                     sparse_rounds: bool = False,
                     output_destination: NewOutputs = None) -> Tuple[list, dict]:
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system.
    There is no issue with using the same interactome interface instance, because they are forked when
    threads are generated and will not interfere.

    :param blank_model_size: the number of uniprots in the blank model
    :param p_val: desired cutoff p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used
     with the number of rounds equal the integer value of that argument
    :param interactome_interface_instance:
    :return: None if no significant nodes, the node and group characteristic
     dictionaries otherwise
    """
    def get_max_for_each_degree(sample_sub_arrray):
        # print('debug max_array_shape:', str(sample_sub_arrray.shape))
        degrees = np.unique(sample_sub_arrray[1, :])
        max_array = []

        for degree in degrees:
            filter = sample_sub_arrray[1, :] == degree
            max_array.append([sample_sub_arrray[0, filter].max(), degree])

        m_arr = np.array(max_array)
        return m_arr.T

    if interactome_interface_instance is None or interactome_interface_instance.node_current == {}:
        raise Exception("tried to compare to blanc an empty interface instance")

    md5_hash = interactome_interface_instance.md5_hash()

    background_sub_array_list = []
    max_sub_array_list = []
    count = 0

    log.info("looking to test against:"
             "\t size: %s \t sys_hash: %s \t sparse_rounds: %s" %
             (blank_model_size, md5_hash, sparse_rounds))

    log.info("samples found to test against:\t %s" %
             count_interactome_rand_samp({'size': blank_model_size,
                                          'sys_hash': md5_hash,
                                          'sparse_rounds': sparse_rounds}))

    background_sample = find_interactome_rand_samp({'size': blank_model_size,
                                                    'sys_hash': md5_hash,
                                                    'sparse_rounds': sparse_rounds})

    for i, sample in enumerate(background_sample):

        _, node_currents = pickle.loads(sample['currents'])

        dict_system = interactome_interface_instance.format_node_props(node_currents, limit=0)
        background_sub_array = list(dict_system.values())

        if np.array(background_sub_array).T.shape[0] < 2:
            log.info(background_sub_array)
            continue

        background_sub_array_list.append(np.array(background_sub_array).T)
        # print(np.array(background_sub_array).T.shape)
        # pprint(background_sub_array)
        max_arr = get_max_for_each_degree(np.array(background_sub_array).T)
        max_sub_array_list.append(max_arr)
        count = i

    # This part declares the pre-operators required for the verification of a
    # real sample

    background_array = np.concatenate(tuple(background_sub_array_list), axis=1)
    max_array = np.concatenate(tuple(max_sub_array_list), axis=1)

    node_currents = interactome_interface_instance.node_current
    dict_system = interactome_interface_instance.format_node_props(node_currents)

    curr_inf_conf_tot = np.array([[int(key)] + list(val) for key, val in list(dict_system.items())]).T

    node_ids, query_array = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1, 2), :])

    log.info("stats on  %s samples" % count)

    background_density = kde_compute(background_array[(1, 0), :], 50, count)
    base_bi_corr = background_array[(0, 1), :]

    r_rels = []
    r_std_nodes = []

    degrees = np.unique(query_array[1, :])

    combined_p_vals = np.ones_like(query_array[1, :])

    for degree in degrees.tolist():
        _filter = query_array[1, :] == degree

        entry = query_array[:, _filter]
        background_set = background_array[:, background_array[1, :] == degree]

        max_current_per_run = get_neighboring_degrees(degree,
                                                      max_array,
                                                      min_nodes=min_nodes_for_p_val)

        p_vals = active_default_p_value_policy(entry, max_current_per_run)
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

    not_random_nodes = [node_id for node_id in node_ids[r_nodes < p_val].tolist()]

    # basically the second element below are the nodes that contribute to the
    #  information flow through the node that is considered as non-random

    log.debug('debug, not random nodes: %s', not_random_nodes)
    log.debug('debug bulbs_id_disp_name: %s',
              list(interactome_interface_instance.neo4j_id_2_display_name.items())[:10])

    node_char_list = [
        [int(nr_node_id), interactome_interface_instance.neo4j_id_2_display_name[nr_node_id]] +
        dict_system[nr_node_id] + r_nodes[node_ids == float(nr_node_id)].tolist()
        for nr_node_id in not_random_nodes]

    nodes_dict = np.hstack((node_ids[:, np.newaxis],
                            r_nodes[:, np.newaxis],
                            r_rels[:, np.newaxis],
                            r_std_nodes[:, np.newaxis]))
    nodes_dict = dict((node[0], (node[1], node[2], node[3])) for node in nodes_dict.tolist())
    nodes_dict = defaultdict(lambda: (1., 0., 0.), nodes_dict)  # corresponds to the cases of super low flow - never significant

    # TODO: pull the groups corresponding to non-random associations.
    # => Will not implement, it's already done by Gephi

    return sorted(node_char_list, key=lambda x: x[4]), nodes_dict


# TODO: [weighted inputs] add support for a dict as source_list, not only list
def auto_analyze(source_list: List[List[int]],
                 output_destinations_list: Union[List[str], None] = None,
                 desired_depth: int = 24,
                 processors: int = 0,
                 background_list: Union[List[int], None] = None,
                 skip_sampling: bool = False,
                 p_value_cutoff: float = -1,
                 ) -> None:
    """
    Automatically analyzes the itneractome synergetic action of the RNA_seq results

    :param source_list: python list of hits for each condition
    :param output_destinations_list: list of names for each condition
    :param desired_depth: total samples we would like to compare each set of hits with
    :param processors: number of processes that will be loaded. as a rule of thumb,
    for max performance, use N-1 processors, where N is the number of physical cores on the
    machine, which is the default
    :param background_list list of physical entities that an experimental method can retrieve
    :param skip_sampling: if true, will skip background sampling step
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

    # TODO: [Better Sampling]
    #  check MongoDb to see if we have enough samples of the needed type, adjust the sampling
    # noinspection PyTypeChecker
    if desired_depth % processors != 0:
        desired_depth = desired_depth // processors + 1
    else:
        desired_depth = desired_depth // processors

    if p_value_cutoff <= 0:
        p_value_cutoff = default_p_val_cutoff

    for hits_list, output_destination in zip(source_list, output_destinations_list):

        log.info('Auto analyzing list of interest: %s', len(hits_list))

        outputs_subdirs = NewOutputs(output_destination)

        interactome_interface = get_interactome_interface(background_up_ids=background_list)

        interactome_interface.set_uniprot_source(list(hits_list))
        log.debug(" e_p_u_b_i length after UP_source was set: %s",
                  len(interactome_interface.active_up_sample))

        if not skip_sampling:
            log.info("spawning a sampler for %s proteins @ %s compops/sec",
                     len(interactome_interface.active_up_sample), estimated_comp_ops)

        # dense analysis
        if len(interactome_interface.active_up_sample) < sparse_analysis_threshold:

            if not skip_sampling:
                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(interactome_interface.active_up_sample),
                         'full',
                         len(interactome_interface.active_up_sample) ** 2 /
                         estimated_comp_ops / 60)

                spawn_sampler_pool(
                    processors,
                    [len(interactome_interface.active_up_sample)],
                    [desired_depth],
                    background_set=background_list)

            interactome_interface.compute_current_and_potentials()

            nr_nodes, p_val_dict = compare_to_blank(
                len(interactome_interface.active_up_sample),
                interactome_interface,
                p_val=p_value_cutoff,
                output_destination=outputs_subdirs
            )

        # sparse analysis
        else:
            ceiling = min(205, len(interactome_interface.active_up_sample))
            sampling_depth = max((ceiling - 5) ** 2 //
                                 len(interactome_interface.active_up_sample),
                                 5)

            if not skip_sampling:

                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(interactome_interface.active_up_sample),
                         sampling_depth,
                         len(interactome_interface.active_up_sample) *
                         sampling_depth / 2 / 60 / estimated_comp_ops)

                spawn_sampler_pool(processors,
                                   [len(interactome_interface.active_up_sample)],
                                   [desired_depth],
                                   sparse_rounds=sampling_depth,
                                   background_set=background_list)

            log.info('real run characteristics: sys_hash: %s, size: %s, sparse_rounds: %s' %
                     (interactome_interface.md5_hash(),
                      len(interactome_interface.active_up_sample), sampling_depth))

            interactome_interface.compute_current_and_potentials(sparse_samples=sampling_depth)

            nr_nodes, p_val_dict = compare_to_blank(
                len(interactome_interface.active_up_sample),
                interactome_interface,
                p_val=p_value_cutoff,
                sparse_rounds=sampling_depth,
                output_destination=outputs_subdirs
            )

        interactome_interface.export_conduction_system(p_val_dict,
                                                       output_location=outputs_subdirs.Interactome_GDF_output)

        # # old results print-out
        # log.info('\t %s \t %s \t %s \t %s \t %s', 'node id',
        #     'display name', 'info flow', 'degree', 'p value')
        #
        # for node in nr_nodes:
        #     log.info('\t %s \t %s \t %.3g \t %d \t %.3g', *node)

        with open(outputs_subdirs.interactome_network_output, 'wt') as output:
            writer = csv_writer(output, delimiter='\t')
            writer.writerow(['node id', 'display name', 'info flow', 'degree', 'p value'])
            for node in nr_nodes:
                writer.writerow(node)

        # using tabulate

        headers = ['node id', 'display name', 'info flow', 'degree', 'p value']

        print(tabulate(nr_nodes, headers, tablefmt='simple', floatfmt=".3g"))


if __name__ == "__main__":

    # pprinter = PrettyPrinter(indent=4)
    # background_set = MatrixGetter(True, False)
    # background_set.fast_load()

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

    local_matrix = InteractomeInterface()
    local_matrix.fast_load()
    # spawn_sampler_pool(3, [50], [3], background_set=None)
    spawn_sampler(([50], [3], False, None, 0))

    # local_matrix.randomly_sample([195], [10], sparse_rounds=195)
