"""
New analytical routines for the interactome
"""
import pickle
from csv import writer as csv_writer
from multiprocessing import Pool
from collections import defaultdict
import traceback
import psutil
from typing import Union, Tuple, List
import numpy as np
from matplotlib import pyplot as plt
from tabulate import tabulate

from bioflow.configs.main_configs import NewOutputs, \
    sparse_analysis_threshold, default_p_val_cutoff, min_nodes_for_p_val, \
    default_background_samples
from bioflow.sample_storage.mongodb import find_interactome_rand_samp, count_interactome_rand_samp
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.general_utils import _is_int
from bioflow.utils.log_behavior import get_logger
from bioflow.neo4j_db.db_io_routines import translate_reweight_dict
from bioflow.algorithms_bank.flow_significance_evaluation import get_neighboring_degrees,\
    get_p_val_by_gumbel
import bioflow.algorithms_bank.sampling_policies as sampling_policies


log = get_logger(__name__)


def get_interactome_interface(background_up_ids=()) -> InteractomeInterface:
    """
    Retrieves an "InteractomeInterface" object

    :return:
    """
    interactome_interface_instance = InteractomeInterface(background_up_ids=background_up_ids)
    interactome_interface_instance.fast_load()
    log.debug("get_interactome state e_p_u_b_i length: %s",
              len(interactome_interface_instance._background))
    # log.info("interactome interface loaded in %s" % interactome_interface_instance.pretty_time())
    # is the case now
    return interactome_interface_instance


# REFACTOR: [MAINTAINABILITY] args puck need to be a named tuple, preferably a typed one
def _spawn_sampler(args_puck):
    """
    Spawns a sampler initialized from the default GO_Interface.

    :param args_puck: combined list of sample sizes, samples to imitate, background sets, and sparse
    sampling argument
    """
    background_set_arg = args_puck[3]

    if args_puck[4] is not None:
        interactome_interface_instance = args_puck[4]
    else:
        interactome_interface_instance = get_interactome_interface(background_set_arg)

    hits_list, sec_list = args_puck[0]
    interactome_interface_instance.set_flow_sources(hits_list, sec_list)

    iterations = args_puck[1]
    sparse_rounds = args_puck[2]

    sampling_policy = args_puck[5]
    sampling_options = args_puck[6]

    pool_no = args_puck[-1]

    interactome_interface_instance.reset_thread_hex()
    interactome_interface_instance.randomly_sample(
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
        forced_interactome_interface,
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
    :param forced_interactome_interface: a provided InteractomeInterface that contains sets to
        imitate
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
                forced_interactome_interface,
                sampling_policy,
                sampling_options)]

    payload_list = payload * pool_size
    payload_list = [list(item) + [i] for i, item in enumerate(payload_list)]  # prepare the payload

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
            pool.terminate()

        # log.info('Pool terminated')

    else:
        log.debug('spawning single-thread sampler with payload %s', payload)
        for _payload in payload_list:
            _spawn_sampler(_payload)


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
    plt.savefig(save_path.interactome_network_scatterplot)
    plt.clf()


def compare_to_blank(interactome_interface_instance: InteractomeInterface,
                     p_val_cutoff: float = 0.05,
                     sparse_rounds: int = -1,
                     output_destination: NewOutputs = None,
                     random_sampling_method=sampling_policies.matched_sampling,
                     random_sampling_option='exact',
                     ) -> Tuple[list, dict]:
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system.
    There is no issue with using the same interactome interface instance, because they are forked when
    threads are generated and will not interfere.

    :param p_val_cutoff: desired cutoff p_value for the returned terms
    :param sparse_rounds: if set to a number, sparse computation technique would be used
        with the number of rounds equal the integer value of that argument
    :param interactome_interface_instance: Interactome interface with loaded real hits list to
        analyse
    :param output_destination: configs object from main_configs, specifying where the results
        will be saved
    :param random_sampling_method: sampling policy used
    :param random_sampling_option: sampling policy optional argument
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
    active_sample_hash = interactome_interface_instance.active_sample_md5_hash(sparse_rounds)

    background_sub_array_list = []
    max_sub_array_list = []
    count = 0

    log.info("looking to test against:"
             "\t target_hash: %s \t sys_hash: %s \n"
             "random sampled according to %s/%s" %
             (active_sample_hash, md5_hash, random_sampling_method.__name__, random_sampling_option))

    samples_to_test_against = count_interactome_rand_samp({
                                          'active_sample_hash': active_sample_hash,
                                          'sys_hash': md5_hash,
                                          'sampling_policy': random_sampling_method.__name__,
                                          'sampling_policy_options': random_sampling_option})

    if samples_to_test_against == 0:
        raise Exception('No samples found to test against. '
                        'There is likely a discrepancy between the parameters used for sampling '
                        'and parameters used for the life sample analysis')

    log.info("samples found to test against:\t %d" % samples_to_test_against)

    background_samples = find_interactome_rand_samp({
                                          'active_sample_hash': active_sample_hash,
                                          'sys_hash': md5_hash,
                                          'sampling_policy': random_sampling_method.__name__,
                                          'sampling_policy_options': random_sampling_option})

    for i, sample in enumerate(background_samples):

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

        # REFACTOR: [MODULARITY] this part is too coupled. we should factor it out
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

    not_random_nodes = [node_id for node_id in node_ids[r_nodes < p_val_cutoff].tolist()]

    # basically the second element below are the nodes that contribute to the
    #  information flow through the node that is considered as non-random

    log.debug('debug, not random nodes: %s', not_random_nodes)
    log.debug('debug bulbs_id_disp_name: %s',
              list(interactome_interface_instance.neo4j_id_2_display_name.items())[:10])

    node_char_list = [
        [int(nr_node_id),
         interactome_interface_instance.neo4j_id_2_display_name[nr_node_id]] +
        dict_system[nr_node_id] + r_nodes[node_ids == float(nr_node_id)].tolist()
        for nr_node_id in not_random_nodes]

    nodes_dict = np.hstack((node_ids[:, np.newaxis],
                            r_nodes[:, np.newaxis],
                            r_rels[:, np.newaxis],
                            r_std_nodes[:, np.newaxis]))
    nodes_dict = dict((node[0], (node[1], node[2], node[3])) for node in nodes_dict.tolist())
    nodes_dict = defaultdict(lambda: (1., 0., 0.), nodes_dict)  # corresponds to the cases of super low flow - never significant

    # NOFX: [structure analysis] pull the groups corresponding to non-random associations.
    # TODO: [total flow] pull out the internode tension map
    # => Will not implement, it's already done by Gephi

    return sorted(node_char_list, key=lambda x: x[4]), nodes_dict


def auto_analyze(source_list: List[Union[List[int], List[Tuple[int, float]]]],
                 secondary_source_list: List[Union[List[int],
                                                   List[Tuple[int, float]],
                                                   None]] = None,
                 output_destinations_list: Union[List[str], None] = None,
                 random_samples_to_test_against: int = default_background_samples,
                 processors: int = 0,
                 background_list: List[Union[List[int], List[Tuple[int, float]], None]] = None,
                 skip_sampling: bool = False,
                 p_value_cutoff: float = -1,
                 sampling_policy=sampling_policies.matched_sampling,
                 sampling_policy_options='exact',
                 explicit_interface=None,
                 forced_lapl_reweight=None,
                 ) -> None:
    """
    Automatically analyzes the interactome synergetic action of the experimental hit lists

    :param source_list: python list of hits for each condition
    :param secondary_source_list: secondary list to which calculate the flow from hist, if needed
    :param output_destinations_list: list of names for each condition
    :param random_samples_to_test_against: total samples we would like to compare each set of hits with
    :param processors: number of processes that will be loaded. as a rule of thumb,
        for max performance, use N-1 processors, where N is the number of physical cores on the
        machine, which is the default
    :param background_list:  list of physical entities that an experimental method can retrieve,
        optionally with weights indicating the likelihood of retrieval at random
    :param skip_sampling: if true, will skip background sparse_sampling step
    :param p_value_cutoff: highest p_value up to which to report the results
    :param sampling_policy: sampling policy used
    :param sampling_policy_options: sampling policy optional argument
    :param explicit_interface: an explicit InteractomeInterface instance in case any of the deep
        defaults (eg flow calculation function) are modified
    :param forced_lapl_reweight: dictionary providing instructions for the modification of
        interactome laplacians weight edges will be set to a given value, nodes will have all the
        edges connecting to them multiplied by the value
    :return:
    """
    # Multiple re-spawns of threaded processing are incompatbile with scikits.sparse.cholmod
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

    if forced_lapl_reweight is not None:
        # print('<<<<<')
        forced_lapl_reweight = translate_reweight_dict(forced_lapl_reweight)
        log.debug("forced reweight instructions dict: %s " % str(forced_lapl_reweight))
        # print('>>>>>')

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

        interactome_interface = get_interactome_interface(background_up_ids=background_list)

        if explicit_interface is not None:
            interactome_interface = explicit_interface

        interactome_interface.set_flow_sources(hits_list, sec_list)

        total_ops = interactome_interface.evaluate_ops()
        sparse_rounds = interactome_interface.reduce_ops(sparse_analysis_threshold**2)

        if sparse_rounds > 0:
            log.info('estimated ops for dense sampling would have been %.1f, '
                     'which is more than the threshold (%d^2). '
                     'Using sparse sampling with %d rounds' % (total_ops,
                                                               sparse_analysis_threshold,
                                                               sparse_rounds))

        else:
            log.info('estimated ops for dense sampling %.1f' % (total_ops))

        md5_hash = interactome_interface.md5_hash()
        active_sample_hash = interactome_interface.active_sample_md5_hash(sparse_rounds)

        in_storage = count_interactome_rand_samp({'active_sample_hash': active_sample_hash,
                                                  'sys_hash': md5_hash,
                                                  'sampling_policy': sampling_policy.__name__,
                                                  'sampling_policy_options': sampling_policy_options})

        if in_storage >= random_samples_to_test_against:
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
                                forced_interactome_interface=explicit_interface,
                                sampling_policy=sampling_policy,
                                sampling_options=sampling_policy_options)

        if forced_lapl_reweight is not None:
            interactome_interface.apply_reweight_dict(forced_lapl_reweight)

        interactome_interface.compute_current_and_potentials()

        nr_nodes, p_val_dict = compare_to_blank(
            interactome_interface,
            p_val_cutoff=p_value_cutoff,
            sparse_rounds=sparse_rounds,
            output_destination=outputs_subdirs,
            random_sampling_method=sampling_policy,
            random_sampling_option=sampling_policy_options
        )

        interactome_interface.export_conduction_system(p_val_dict,
                                                       output_location=outputs_subdirs.Interactome_GDF_output)

        with open(outputs_subdirs.interactome_network_output, 'wt') as output:
            writer = csv_writer(output, delimiter='\t')
            writer.writerow(['node id', 'display name', 'info flow', 'degree', 'p value'])
            for node in nr_nodes:
                writer.writerow(node)

        # using tabulate

        headers = ['node id', 'display name', 'info flow', 'degree', 'p value']

        print(tabulate(nr_nodes, headers, tablefmt='simple', floatfmt=".3g"))


if __name__ == "__main__":

    local_matrix = InteractomeInterface()
    local_matrix.fast_load()