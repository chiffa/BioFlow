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

import numpy as np
from matplotlib import pyplot as plt

from bioflow.algorithms_bank.conduction_routines import perform_clustering
from bioflow.main_configs import interactome_rand_samp_db, Outputs, Dumps, estimated_comp_ops
from bioflow.molecular_network.InteractomeInterface import InteractomeInterface
from bioflow.utils.dataviz import kde_compute
from bioflow.utils.log_behavior import get_logger
from bioflow.utils.io_routines import get_source_bulbs_ids, get_background_bulbs_ids

from scipy.stats import gumbel_r

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
              len(interactome_interface_instance.entry_point_uniprots_neo4j_ids))
    log.info("interactome interface loaded in %s" % interactome_interface_instance.pretty_time())
    # is the case now
    return interactome_interface_instance


# TODO: defined only as an auxilary to the spawn_sampler_pool => inline that in the sampler pool
def spawn_sampler(args_puck):
    """
    Spawns a sampler initialized from the default GO_Interface.

    :param args_puck: combined list of sample sizes and
    iterations (required for Pool.map usage)
    """
    interactome_interface_instance_arg = args_puck[3]

    if interactome_interface_instance_arg is None:
        interactome_interface_instance = get_interactome_interface()
    else:
        interactome_interface_instance = interactome_interface_instance_arg

    sample_size_list = args_puck[0]
    iteration_list = args_puck[1]
    sparse_rounds = args_puck[2]
    pool_no = args_puck[-1]

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
        interactome_interface_instance,
        sparse_rounds=False):
    """
    Spawns a pool of samplers of the information flow within the GO system

    :param pool_size: number of processes that are performing the sample pooling and analyzing
    :param sample_size_list: size of the sample list
    :param interaction_list_per_pool: number of iterations performing the pooling of the samples
     in each list
    :param sparse_rounds:
    :param interactome_interface_instance:
    """
    try:
        process_pool = Pool(pool_size)
        payload = [
            (sample_size_list,
             interaction_list_per_pool,
             sparse_rounds,
             interactome_interface_instance)]
        log.debug('spawning the sampler with payload %s', payload)
        payload_list = payload * pool_size
        payload_list = [list(item)+[i] for i, item in enumerate(payload_list)]
        process_pool.map(spawn_sampler, payload_list)

    except Exception as e:
        msg = "{}\n\nOriginal {}".format(e, traceback.format_exc())
        raise type(e)(msg)


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


def samples_scatter_and_hist(background_curr_deg_conf, true_sample_bi_corr_array):
    """
    A general function that performs demonstration of an example of random samples of
     the same size as our sample and of our sample and conducts the statistical tests
     on wherther any of nodes or functional groups in our sample are non-random

    :param background_curr_deg_conf: [[current, informativity, confusion_potential], ...] -
    characteristics of the random samples
    :param true_sample_bi_corr_array: [[current, informativity, confusion_potential], ...] -
    characteristics of the true sample. If none, nothing happens
    :return: None
    """

    fig = plt.figure()
    fig.set_size_inches(30, 20)

    plt.subplot(211)
    plt.title('current through nodes')

    bins = np.linspace(
        background_curr_deg_conf[0, :].min(),
        background_curr_deg_conf[0, :].max(), 100)

    if true_sample_bi_corr_array is not None:
        bins = np.linspace(min(background_curr_deg_conf[0, :].min(), true_sample_bi_corr_array[0, :].min()),
                           max(background_curr_deg_conf[0, :].max(), true_sample_bi_corr_array[0, :].max()),
                           100)

    plt.hist(background_curr_deg_conf[0, :],
             bins=bins, histtype='step', log=True, color='b')

    if true_sample_bi_corr_array is not None:
        plt.hist(true_sample_bi_corr_array[0, :],
                 bins=bins, histtype='step', log=True, color='r')

    plt.subplot(212)
    plt.scatter(background_curr_deg_conf[1, :], background_curr_deg_conf[0, :])

    if true_sample_bi_corr_array is not None:
        plt.scatter(true_sample_bi_corr_array[1, :], true_sample_bi_corr_array[0, :],
                    color='r', alpha=0.5)

    # plt.show()
    plt.savefig(Outputs.interactome_network_stats)
    plt.clf()



def compare_to_blank(blank_model_size, p_val=0.05, sparse_rounds=False,
                     interactome_interface_instance=None):
    """
    Recovers the statistics on the circulation nodes and shows the visual of a circulation system.
    There is no issue with using the same interactome interface instance, because they are forked when
    threads are generated and will not interfere.

    :param blank_model_size: the number of uniprots in the blank model
    :param p_val: desired p_value for the returned terms
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

    if interactome_interface_instance is None:
        interactome_interface_instance = InteractomeInterface(True, True)
        interactome_interface_instance.fast_load()

    md5_hash = interactome_interface_instance.md5_hash()

    background_sub_array_list = []
    max_sub_array_list = []
    count = 0

    log.info("looking to test against:"
             "\t size: %s \t sys_hash: %s \t sparse_rounds: %s" %
             (blank_model_size, md5_hash, sparse_rounds))

    log.info("samples found to test against:\t %s" %
             interactome_rand_samp_db.find({'size': blank_model_size,
                                            'sys_hash': md5_hash,
                                            'sparse_rounds': sparse_rounds}).count())

    for i, sample in enumerate(interactome_rand_samp_db.find(
                                                            {'size': blank_model_size,
                                                             'sys_hash': md5_hash,
                                                             'sparse_rounds': sparse_rounds})):

        _, node_currents = pickle.loads(sample['currents'])

        dictionary_system = interactome_interface_instance.format_node_props(node_currents, limit=0)
        background_sub_array = list(dictionary_system.node_current_values())
        if np.array(background_sub_array).T.shape[0] < 2:
            print(background_sub_array)
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
    dictionary_system = interactome_interface_instance.format_node_props(node_currents)
    curr_inf_conf_tot = np.array(
        [[int(key)] + list(val) for key, val in list(dictionary_system.items())]).T

    node_ids, query_array = (curr_inf_conf_tot[0, :], curr_inf_conf_tot[(1, 2), :])

    log.info("stats on  %s samples" % count)

    background_density = kde_compute(background_array[(1, 0), :], 50, count)
    base_bi_corr = background_array[(0, 1), :]

    r_rels = []
    r_std_nodes = []

    # TODO: idea for the improved statistics, cluster a test node of degree k with 100 nodes with
    #  closest degrees

    samples_scatter_and_hist(background_array, query_array)

    degrees = np.unique(query_array[1, :])

    combined_p_vals = np.ones_like(query_array[1, :])

    for degree in degrees.tolist():
        filter = query_array[1, :] == degree

        entry = query_array[:, filter]
        background_set = background_array[:, background_array[1, :] == degree]
        max_set = max_array[:, max_array[1, :] == degree]

        params = gumbel_r.fit(max_set[0, :])
        arg = params[:-2]
        mu = params[-2]
        beta = params[-1]

        frozen_gumbel = gumbel_r(loc=mu, scale=beta)

        p_vals = 1 - frozen_gumbel.cdf(entry[0, :])

        combined_p_vals[filter] = p_vals

        # TODO: insert into appropriate locations => we assume that the order is preserved

        # samples_scatter_and_hist(max_set, entry)

    r_nodes = background_density(query_array[(1, 0), :])  # this is currently used as a p-value, which is problematic.
    r_nodes = combined_p_vals

    for point in query_array.T:
        selector = np.logical_and(base_bi_corr[1, :] > point[1]*0.9, base_bi_corr[1, :] < point[1]*1.1)
        r_rels.append(point[0]/np.mean(base_bi_corr[0, selector]))
        r_std_nodes.append((point[0]-np.mean(base_bi_corr[0, selector]))/np.std(base_bi_corr[0, selector]))

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
        dictionary_system[nr_node_id] + r_nodes[node_ids == float(nr_node_id)].tolist()
        for nr_node_id in not_random_nodes]

    nodes_dict = np.hstack((node_ids[:, np.newaxis], r_nodes[:, np.newaxis], r_rels[:, np.newaxis], r_std_nodes[:, np.newaxis]))
    nodes_dict = dict((node[0], (node[1], node[2], node[3])) for node in nodes_dict.tolist())
    nodes_dict = defaultdict(lambda: (1., 0., 0.), nodes_dict)  # corresponds to the cases of super low flow - never significant

    # TODO: pull the groups corresponding to non-random associations.

    return sorted(node_char_list, key=lambda x: x[4]), nodes_dict



# TODO: source_list is a single list, not list of lists.
# TODO: add support for a dict as source_list, not only list
def auto_analyze(source_list,
                 desired_depth=24, processors=4,
                 background_list=None, skip_sampling=False,
                 from_memoization=False):
    """
    Automatically analyzes the itneractome synergetic action of the RNA_seq results

    :param source_list: python list of "hit" physical entities
    :param desired_depth: desired sampling depth
    :param processors: number of processes that will be loaded. as a rule of thumb,
    for max performance, use N-1 processors, where N is the number of physical cores on the machine
    :param background_list list of physical entities that an experimental method can retrieve
    :param skip_sampling: if true, will skip background sampling step
    """
    # noinspection PyTypeChecker
    if desired_depth % processors != 0:
        desired_depth = desired_depth // processors + 1
    else:
        desired_depth = desired_depth // processors

    for _list in source_list:
        log.info('Auto analyzing list of interest: %s', len(_list))

        interactome_interface = get_interactome_interface()
        log.debug("retrieved interactome_interface instance e_p_u_b_i length: %s",
                  len(interactome_interface.entry_point_uniprots_neo4j_ids))

        interactome_interface.set_uniprot_source(list(_list))
        log.debug(" e_p_u_b_i length after UP_source was set: %s",
                  len(interactome_interface.entry_point_uniprots_neo4j_ids))

        if background_list:
            interactome_interface.background = background_list
            interactome_interface.connected_uniprots = list(
                set(interactome_interface.connected_uniprots).intersection(set(interactome_interface.background)))

        if not skip_sampling:
            log.info("spawning a sampler for %s proteins @ %s compops/sec",
                     len(interactome_interface.entry_point_uniprots_neo4j_ids), estimated_comp_ops)

        if len(interactome_interface.entry_point_uniprots_neo4j_ids) < 60:

            if not skip_sampling:
                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(interactome_interface.entry_point_uniprots_neo4j_ids),
                         'full',
                         len(interactome_interface.entry_point_uniprots_neo4j_ids) ** 2 /
                         estimated_comp_ops / 60)

                spawn_sampler_pool(
                    processors,
                    [len(interactome_interface.entry_point_uniprots_neo4j_ids)],
                    [desired_depth],
                    interactome_interface_instance=interactome_interface)

            if not from_memoization:
                interactome_interface.build_extended_conduction_system()
            else:
                interactome_interface.build_extended_conduction_system(fast_load=True)

            nr_nodes, p_val_dict = compare_to_blank(
                len(interactome_interface.entry_point_uniprots_neo4j_ids),
                p_val=0.9,
                interactome_interface_instance=interactome_interface)


        else:
            ceiling = min(205, len(interactome_interface.entry_point_uniprots_neo4j_ids))
            sampling_depth = max((ceiling-5) ** 2 /
                                 len(interactome_interface.entry_point_uniprots_neo4j_ids),
                                 5)

            if not skip_sampling:
                log.info('length: %s \t sampling depth: %s \t, estimated round time: %s min',
                         len(interactome_interface.entry_point_uniprots_neo4j_ids),
                         sampling_depth,
                         len(interactome_interface.entry_point_uniprots_neo4j_ids) *
                         sampling_depth / 2 / 60)

                spawn_sampler_pool(processors,
                                   [len(interactome_interface.entry_point_uniprots_neo4j_ids)],
                                   [desired_depth],
                                   sparse_rounds=sampling_depth,
                                   interactome_interface_instance=interactome_interface)

            log.info('real run characteristics: sys_hash: %s, size: %s, sparse_rounds: %s' % (interactome_interface.md5_hash(),
                                                                                              len(interactome_interface.entry_point_uniprots_neo4j_ids), sampling_depth))

            if not from_memoization:
                interactome_interface.build_extended_conduction_system(sparse_samples=sampling_depth)
                # interactome_interface.export_conduction_system()
            else:
                interactome_interface.build_extended_conduction_system(sparse_samples=sampling_depth, fast_load=True)
            nr_nodes, p_val_dict = compare_to_blank(
                len(interactome_interface.entry_point_uniprots_neo4j_ids), p_val=0.9,
                sparse_rounds=sampling_depth, interactome_interface_instance=interactome_interface)

        interactome_interface.export_conduction_system(p_val_dict)

        log.info('\t %s \t %s \t %s \t %s \t %s', 'node id',
            'display name', 'info flow', 'degree', 'p value')

        for node in nr_nodes:
            log.info('\t %s \t %s \t %s \t %s \t %s', *node)

        with open(Outputs.interactome_network_output, 'wt') as output:
            writer = csv_writer(output, delimiter='\t')
            writer.writerow(['node id', 'display name', 'info flow', 'degree', 'p value'])
            for node in nr_nodes:
                writer.writerow(node)


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

    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
    local_matrix.fast_load()

    # spawn_sampler_pool(3, [2], [150], interactome_interface_instance=None)

    # local_matrix.randomly_sample([195], [10], sparse_rounds=195)
